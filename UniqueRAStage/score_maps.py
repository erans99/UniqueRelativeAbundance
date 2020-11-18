import glob
import math
import os
import pickle
import time
from functools import reduce
import logging
import numpy
import pandas

_log = logging.getLogger('URA')

MIN_PARTS_TO_ESTIMATE = 5


def get_stats_single(binsfile, st_perfect):
    bins_dict = pickle.load(open(binsfile, 'rb'), encoding='bytes')
    stats_perfect = {}
    num_parts = 0
    contigs = bins_dict.keys()
    for c in contigs:
        c_perfect = st_perfect[st_perfect.contigs == c]
        bins = bins_dict[c]
        n_bins = (len(bins) - 1)
        num_parts += n_bins
        if len(c_perfect) == 0:
            stats_perfect[int(c)] = [0] * n_bins
            continue
        poses = map(lambda x: int(x), c_perfect.poses.values)
        bin_poses = numpy.digitize(list(poses), bins, right=True)
        vals, cnts = numpy.unique(bin_poses, return_counts=True)
        bin_fills = numpy.zeros(n_bins+2)
        numpy.put(bin_fills, vals, cnts)
        stats_perfect[int(c)] = list(bin_fills[1:-1])
    return stats_perfect, num_parts, len(st_perfect)


def get_stats(filename, bins_dict_files):
    all_unique = pandas.read_pickle(filename)
    _log.info("Starting with %d reads" % len(all_unique))
    strains = all_unique.IDs.value_counts().index
    stats_perfect = {}
    num_parts = {}
    num_perfect = {}
    for st in strains:
        binsfile = glob.glob(bins_dict_files % (st + '.pkl'))
        if len(binsfile) == 1:
            binsfile = binsfile[0]
        else:
            _log.warning("Failed reading info file %s %s (%d)" % (st, os.path.join(bins_dict_files % (st + '.pkl')),
                                                           len(binsfile)))
            raise Exception("Invalid bins file %s %s (%d)" % (st, os.path.join(bins_dict_files % (st + '.pkl')),
                                                              len(binsfile)))
        stats_perfect[st], num_parts[st], num_perfect[st] = get_stats_single(binsfile,
                                                                             all_unique[all_unique.IDs == st])
    return strains, stats_perfect, num_perfect, num_parts


# search for a good enough number of bins, by "lion in the dessert" type search
def calc_num_bins(bins_dict, old_needed, exp_th, max_steps=20):
    full_fills = reduce(lambda x, y: x[:-1] + y, bins_dict.values())[:-1]
    res = pandas.Series()
    fills = numpy.array(full_fills)
    # shuffle break potential plasmids/horizontal transfers apart. Is that what we want?
    #    random.shuffle(fills)
    needed_parts = old_needed
    rng = [-1, -1]
    while True:
        if (needed_parts in res.index) or (len(res) > max_steps):
            best = res[res > exp_th]
            if len(best) > 0:
                best = best.idxmin()
                if res.loc[best] < 2 * exp_th:
                    # print("Need %d parts to overpass threshold (%d tries)" % (best, len(res) + 1))
                    return best
            best = res[res <= exp_th]
            if len(best) > 0:
                best = best.idxmax()
                if res.loc[best] > 0.5 * exp_th:
                    # print("Need %d parts to underpass threshold (%d tries)" % (best, len(res) + 1))
                    return best
            # print("Can't find number of needed parts (%d tries)" % len(res))
            return 0
        take_len = int(len(fills) / needed_parts)
        if take_len == 0:
            res.loc[needed_parts] = -1
            needed_parts = len(fills)
            continue
        accum_fills = fills[:take_len * needed_parts].reshape([take_len, needed_parts]).sum(1)
        accum_fills.sort()
        median = accum_fills[int(take_len / 2)]
        if (median > exp_th) and ((median < 1.5 * exp_th) or (needed_parts == 1)):
            # print("Need %d parts to pass threshold (%d tries)" % (needed_parts, len(res) + 1))
            return needed_parts
        res.loc[needed_parts] = median
        if median < exp_th:
            rng[0] = max(rng[0], needed_parts)
            if rng[1] == -1:
                needed_parts *= 2
            else:
                needed_parts = int((needed_parts + rng[1]) / 2)
        else:
            if rng[1] == -1:
                rng[1] = needed_parts
            else:
                rng[1] = min(rng[1], needed_parts)
            if rng[0] == -1:
                needed_parts = max(int(needed_parts / 2), 1)
            else:
                needed_parts = int((needed_parts + rng[0]) / 2)


def calc_needed_parts(stats_perfect, exp_mean, exp_th):
    if exp_mean == 0:
        return 0
    if exp_mean < exp_th:
        old_needed = int(math.ceil((1.0 * exp_th) / exp_mean))
    else:
        old_needed = 1
    needed_parts = calc_num_bins(stats_perfect, old_needed, exp_th)
    # print("Got %d needed parts, vs %d in old method" % (needed_parts, old_needed))
    return needed_parts


def analyse_stats_single(needed_parts, min_partial, stats_perfect, scoresfile, num_uniq):
    num_uniq = num_uniq * 1.
    scores_dict = pickle.load(open(scoresfile, 'rb'), encoding='bytes')
    min_exp = int(min_partial * needed_parts * num_uniq)
    norm_all_fill = []
    num_not_used = 0
    for c in scores_dict.keys():
        fill = 0
        expect = 0
        cnt_cs = 0
        for b, fill_b in enumerate(stats_perfect[c]):
            fill += fill_b
            expect += scores_dict[c][b]
            cnt_cs += 1
            if cnt_cs == needed_parts:
                if expect > min_exp:
                    norm_all_fill.append((num_uniq * float(fill)) / expect)
                elif expect != 0:
                    num_not_used += cnt_cs
                fill = 0
                expect = 0
                cnt_cs = 0
        if cnt_cs > 0:
            if expect > min_exp:
                norm_all_fill.append((num_uniq * float(fill)) / expect)
            elif expect != 0:
                num_not_used += cnt_cs

    if len(norm_all_fill) == 0:
        # print("  on blocks of %d parts (0 parts, %d bins not used)" % (needed_parts, num_not_used))
        return None
    if len(norm_all_fill) < MIN_PARTS_TO_ESTIMATE:
        # s = stats.describe(norm_all_fill)
        # print("  on blocks of %d parts (%d parts, %d bins not used): mean %g (can't do dense mean or variance)" % \
        #       (needed_parts, len(norm_all_fill), num_not_used, s.mean))
        return None
    s = dense_non0_mean(norm_all_fill, needed_parts)
    # if s[1] == 0:
    #     print("  on blocks of %d parts (%d parts, %d bins not used): dense starts at 0 can't calculate" %
    #           (needed_parts, len(norm_all_fill), num_not_used))
    # else:
    #     print("  on blocks of %d parts (%d parts, %d bins not used): dense mean %g std ~%g" %
    #           (needed_parts, len(norm_all_fill), num_not_used, s[1], s[2]))
    return [needed_parts, len(norm_all_fill), s[1], s[2]], s[1]


def analyse_stats(outputf, strains, stats_perfect, num_perfect, num_parts, scores_dict_files, exp_th, min_partial,
                  min_abund, num_uniq):
    exp_th = numpy.float(exp_th)
    min_partial = numpy.float(min_partial)
    abnds = {}
    cover = {}
    _log.info("Got something in %d strains" % len(stats_perfect))
    for st in strains:
        scoresfile = glob.glob(scores_dict_files % (st + '.pkl'))
        if len(scoresfile) == 1:
            scoresfile = scoresfile[0]
        else:
            _log.warning("Failed reading info file %s %s (%d)" % (st, os.path.join(scores_dict_files % (st + '.pkl')),
                                                           len(scoresfile)))
            raise Exception("Invalid score file %s %s (%d)" % (st, os.path.join(scores_dict_files % (st + '.pkl')),
                                                               len(scoresfile)))
        needed_parts = calc_needed_parts(stats_perfect[st], ((1.0 * num_perfect[st]) / num_parts[st]), exp_th)
        if needed_parts == 0:
            continue
        res = analyse_stats_single(needed_parts, min_partial, stats_perfect[st], scoresfile, num_uniq)
        if res is not None:
            cover[st], abnds[st] = res

    abnds = pandas.Series(abnds)
    if sum(abnds) > 0:
        abnds /= sum(abnds)
    else:
        return None
    abnds.sort_values(inplace=True, ascending=False)
    abnds[abnds < min_abund] = 0
    abnds /= sum(abnds)
    abnds = abnds[abnds > 0]
    abnds.to_csv(outputf, header=False, compression='gzip')

    cover = pandas.DataFrame(cover, index=['needed_parts', 'num_parts', 'mean_part_cover', 'std_part_cover']).T
    cover.loc[abnds.index].to_csv(outputf.replace("abundance_", "cover_"), header=True, compression='gzip')
    return abnds


def dense_non0_mean(x, parts_binned, max_0s=0.5, min_parts=10, noise_level=0.001):
    # there is a problem here if there is more then one such peak (i.e not normal + few del/copy parts)
    # problem is worse the smaller is q. I wouldn't go under 0.5
    q = 0.5
    fq = 2.65  # std factor for taking dense 50%

    num_0 = x.count(0)
    if (num_0 > max_0s * len(x)) or ((len(x) - num_0) < min_parts):
        return [0, 0], 0, -1
    x1 = numpy.array(x)
    x1 = x1[x1.nonzero()[0]]
    min_non0 = x1.min()
    num_min = sum(x1 == min_non0)
    # in order to avoid equalities
    x1 += numpy.random.normal(0, noise_level / parts_binned, len(x1))
    x1.sort()
    d = int(q * len(x))
    diff = (x1[d:] - x1[:-d])
    if len(diff) == 0:
        # print("Got empty sequence?!? %d %d %d %d" % (len(x), len(x1), d, len(diff)))
        return [0, x1[0]], 0, -1
    min_pos = diff.argmin()
    if min_pos > num_min:
        return [x1[min_pos], x1[min_pos + d]], x1[min_pos:min_pos + d].mean(), fq * (
                x1[min_pos:min_pos + d].var() ** 0.5)
    else:
        return [0, x1[min_pos + d]], 0, -1


def score_map_file(outputf, inputf, bins_dict_files, scores_dict_files, exp_th, min_partial, min_abund, num_uniq):
    _log.info(("File %s started at " % os.path.basename(inputf)) + time.ctime())
    strains, stats_perfect, num_perfect, num_parts = get_stats(inputf, bins_dict_files)
    _log.info("Mid at " + time.ctime())
    analyse_stats(outputf, strains, stats_perfect, num_perfect, num_parts, scores_dict_files, exp_th, min_partial,
                  min_abund, num_uniq)
    _log.info("Ended at " + time.ctime())


def run_sample(sample, out_file, in_path, scores_dict_files, exp_th, min_partial, bins_dict_files, min_abund, num_uniq):
    infile = glob.glob(os.path.join(in_path, sample + '.pkl'))[0]
    score_map_file(out_file, infile, bins_dict_files, scores_dict_files, exp_th, min_partial, min_abund, num_uniq)
    _log.info('Done calculating abundance for sample %s' % sample)
    return out_file
