#ifndef ClusterTimeMerge_H
#define ClusterTimeMerge_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

// Merging of clusters that are close in time, used by ClusterFinder (TimeMergeWindow).
//
// Clusters are ordered by the median time of their hits and chained: a cluster joins the
// current group when its median is less than `window` ns after the median of the cluster
// before it. Clusters whose median is earlier than `min_time` are never merged. The hits of a
// group are combined, skipping hits (same tube and time) already taken from an earlier cluster
// of the group; hits repeated inside one cluster are kept as ClusterFinder made them. The group
// is keyed by the mean time of its hits, the same time convention ClusterFinder uses.
namespace ClusterTimeMerge {

inline double Median(std::vector<double> t) {
  if (t.empty()) return 0.;
  std::sort(t.begin(), t.end());
  size_t n = t.size();
  return (n % 2) ? t[n / 2] : 0.5 * (t[n / 2 - 1] + t[n / 2]);
}

template <typename HitT>
double MedianTime(const std::vector<HitT>& hits) {
  std::vector<double> t;
  t.reserve(hits.size());
  for (const auto& h : hits) t.push_back(h.GetTime());
  return Median(t);
}

// Cluster keys grouped for merging, groups and members in time order; singletons included.
template <typename HitT>
std::vector<std::vector<double>> Groups(const std::map<double, std::vector<HitT>>& clusters,
                                        double window, double min_time) {
  std::vector<std::pair<double, double>> order;  // (median, key)
  for (const auto& c : clusters) order.emplace_back(MedianTime(c.second), c.first);
  std::stable_sort(order.begin(), order.end(),
                   [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                     return a.first < b.first;
                   });
  std::vector<std::vector<double>> groups;
  bool have_prev = false;
  double prev = 0.;
  for (const auto& o : order) {
    if (o.first < min_time) {
      groups.push_back({o.second});
      continue;
    }
    if (!have_prev || o.first - prev >= window) groups.push_back({});
    groups.back().push_back(o.second);
    prev = o.first;
    have_prev = true;
  }
  return groups;
}

// Merges in place; returns the number of clusters removed.
template <typename HitT>
int Merge(std::map<double, std::vector<HitT>>& clusters,
          std::map<double, std::vector<unsigned long>>& detkeys,
          double window, double min_time) {
  int removed = 0;
  for (const auto& group : Groups(clusters, window, min_time)) {
    if (group.size() < 2) continue;
    std::vector<HitT> hits;
    std::vector<unsigned long> keys;
    std::set<std::pair<int, double>> seen;
    for (double key : group) {
      const std::vector<HitT>& chits = clusters.at(key);
      const std::vector<unsigned long>& ckeys = detkeys.at(key);
      std::set<std::pair<int, double>> mine;
      for (size_t i = 0; i < chits.size(); ++i) {
        std::pair<int, double> id(chits[i].GetTubeId(), chits[i].GetTime());
        if (seen.count(id)) continue;
        mine.insert(id);
        hits.push_back(chits[i]);
        keys.push_back(ckeys.at(i));
      }
      seen.insert(mine.begin(), mine.end());
      clusters.erase(key);
      detkeys.erase(key);
    }
    double newkey = 0.;
    for (const auto& h : hits) newkey += h.GetTime();
    newkey /= hits.size();
    while (clusters.count(newkey)) newkey = std::nextafter(newkey, std::numeric_limits<double>::max());
    clusters.emplace(newkey, hits);
    detkeys.emplace(newkey, keys);
    removed += group.size() - 1;
  }
  return removed;
}

}  // namespace ClusterTimeMerge

#endif
