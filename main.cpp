#include <bits/stdc++.h>
#define REP(i, n) for (int i = 0; (i) < (int)(n); ++ (i))
#define REP3(i, m, n) for (int i = (m); (i) < (int)(n); ++ (i))
#define REP_R(i, n) for (int i = (int)(n) - 1; (i) >= 0; -- (i))
#define REP3R(i, m, n) for (int i = (int)(n) - 1; (i) >= (int)(m); -- (i))
#define ALL(x) std::begin(x), std::end(x)
using namespace std;

template <class T> using reversed_priority_queue = priority_queue<T, vector<T>, greater<T> >;

class xor_shift_128 {
public:
    typedef uint32_t result_type;
    xor_shift_128(uint32_t seed = 42) {
        set_seed(seed);
    }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator() () {
        uint32_t t = (a ^ (a << 11));
        a = b; b = c; c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }
    static constexpr uint32_t max() { return numeric_limits<result_type>::max(); }
    static constexpr uint32_t min() { return numeric_limits<result_type>::min(); }
private:
    uint32_t a, b, c, d;
};

constexpr array<int, 4> DIR_Y = {-1, 1, 0, 0};
constexpr array<int, 4> DIR_X = {0, 0, 1, -1};

string convert_to_command_string(const vector<pair<int, int>>& result) {
    assert (not result.empty());
    string ans;
    REP (i, (int)result.size() - 1) {
        auto [ay, ax] = result[i];
        auto [by, bx] = result[i + 1];
        while (true) {
            if (by < ay and bx == ax) {
                ans.push_back('U');
                -- ay;
            } else if (by > ay and bx == ax) {
                ans.push_back('D');
                ++ ay;
            } else if (by == ay and bx > ax) {
                ans.push_back('R');
                ++ ax;
            } else if (by == ay and bx < ax) {
                ans.push_back('L');
                -- ax;
            } else if (by == ay and bx == ax) {
                break;
            } else {
                assert (false);
            }
        }
    }
    return ans;
}

vector<vector<int>> run_dijkstra(int n, int sy, int sx, const vector<vector<int>> &c) {
    vector<vector<int>> dist(n, vector<int>(n, INT_MAX));
    reversed_priority_queue<tuple<int, int, int>> que;
    dist[sy][sx] = 0;
    que.emplace(dist[sy][sx], sy, sx);
    while (not que.empty()) {
        auto [dist_y_x, y, x] = que.top();
        que.pop();
        if (dist[y][x] < dist_y_x) continue;
        REP (dir, 4) {
            int ny = y + DIR_Y[dir];
            int nx = x + DIR_X[dir];
            if (0 <= ny and ny < n and 0 <= nx and nx < n) {
                if (c[ny][nx] != INT_MAX) {
                    if (dist[y][x] + c[ny][nx] < dist[ny][nx]) {
                        dist[ny][nx] = dist[y][x] + c[ny][nx];
                        que.emplace(dist[ny][nx], ny, nx);
                    }
                }
            }
        }
    }
    return dist;
}

template <class T, class RandomEngine>
T get_sample(const vector<T> &xs, RandomEngine &gen) {
    int n = xs.size();
    assert (n >= 1);
    int i = uniform_int_distribution<int>(0, n - 1)(gen);
    return xs[i];
}

template <class RandomEngine>
string solve(const int n, const int start_y, const int start_x, const vector<vector<int>> &original_c, RandomEngine& gen, chrono::high_resolution_clock::time_point clock_end) {
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();

    // Make a copy
    vector<vector<int>> c = original_c;

#ifdef VISUALIZE
    cerr << "input:" << endl;
    REP (y, n) {
        REP (x, n) {
            cerr << (char)(c[y][x] == INT_MAX ? ' ' : '0' + c[y][x]);
        }
        cerr << endl;
    }
#endif

    // Run Dijkstra from (start_y, start_x)
    vector<vector<int>> dist_from_start = run_dijkstra(n, start_y, start_x, c);
    vector<vector<int>> dist_to_start = dist_from_start;
    REP (y, n) {
        REP (x, n) {
            if (c[y][x] != INT_MAX) {
                dist_to_start[y][x] -= c[y][x];
                dist_to_start[y][x] += c[start_y][start_x];
            }
        }
    }

    auto list_neighbors = [&](int y, int x) -> vector<pair<int, int>> {
        vector<pair<int, int>> neighbors;
        REP (dir, 4) {
            int ny = y + DIR_Y[dir];
            int nx = x + DIR_X[dir];
            if (0 <= ny and ny < n and 0 <= nx and nx < n) {
                if (c[ny][nx] != INT_MAX) {
                    neighbors.emplace_back(ny, nx);
                }
            }
        }
        return neighbors;
    };

    // Collect intersections of roads
    int intersection_count = 0;
    vector<vector<int>> intersection_index(n, vector<int>(n, -1));
    REP (y, n) {
        REP (x, n) {
            if (c[y][x] != INT_MAX) {
                auto neighbors = list_neighbors(y, x);
                int sum_dy = 0;
                int sum_dx = 0;
                for (auto [ny, nx] : neighbors) {
                    sum_dy += ny - y;
                    sum_dx += nx - x;
                }
                if (neighbors.size() >= 3 or (neighbors.size() == 2 and (sum_dy != 0 or sum_dx != 0))) {
                    intersection_index[y][x] = intersection_count;
                    intersection_count += 1;
                }
            }
        }
    }
#ifdef LOCAL
    cerr << "intersection count = " << intersection_count << endl;
#endif

    // Make the succinct graph of intersections
    vector<vector<pair<int, int>>> g(intersection_count);
    vector<pair<int, int>> intersection_location(intersection_count);
    REP (y, n) {
        REP (x, n) {
            int i = intersection_index[y][x];
            if (i == -1) {
                continue;
            }
            intersection_location[i] = make_pair(y, x);

            // up and left
            for (int dir : {0, 3}) {
                int ny = y + DIR_Y[dir];
                int nx = x + DIR_X[dir];
                if (not (0 <= ny and ny < n and 0 <= nx and nx < n) or c[ny][nx] == INT_MAX) {
                    continue;
                }
                int cost = 0;
                while (0 <= ny and ny < n and 0 <= nx and nx < n and c[ny][nx] != INT_MAX and intersection_index[ny][nx] == -1) {
                    cost += c[ny][nx];
                    ny += DIR_Y[dir];
                    nx += DIR_X[dir];
                }
                if (not (0 <= ny and ny < n and 0 <= nx and nx < n) or c[ny][nx] == INT_MAX) {
                    continue;
                }
                int j = intersection_index[ny][nx];
                g[i].emplace_back(j, cost + c[ny][nx]);
                g[j].emplace_back(i, cost + c[y][x]);
            }
        }
    }
#ifdef VISUALIZE
    REP (y, n) {
        REP (x, n) {
            int i = intersection_index[y][x];
            if (i == -1) {
                continue;
            }
            for (auto [j, cost] : g[i]) {
                auto [ny, nx] = intersection_location[j];
                // cerr << "(" << y << ", " << x << ") -> (" << ny << ", " << nx << "): " << cost << endl;
            }
        }
    }
#endif

    auto dist_to_start_intersection = [&](int i) {
        assert (0 <= i and i < intersection_count);
        auto [y, x] = intersection_location[i];
        return dist_to_start[y][x];
    };
    auto dist_from_start_intersection = [&](int i) {
        assert (0 <= i and i < intersection_count);
        auto [y, x] = intersection_location[i];
        return dist_from_start[y][x];
    };

    // Warshall-Floyd on the intersections graph
    vector<vector<int>> dist(intersection_count, vector<int>(intersection_count, 1e9));  // use small INF
    REP (i, intersection_count) {
        for (auto [j, cost] : g[i]) {
            dist[i][j] = cost;
        }
    }
    REP (k, intersection_count) {
        REP (i, intersection_count) {
            REP (j, intersection_count) {
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
            }
        }
    }

    // List streets
    int street_count = 0;
    vector<array<int, 2>> streets_of(intersection_count, {-1, -1});
    REP (y, n) {
        REP (x, n) {
            // up and left
            for (int dir : {0, 3}) {
                int ny = y + DIR_Y[dir];
                int nx = x + DIR_X[dir];
                // the upper end or the left end of a street
                if (not (0 <= ny and ny < n and 0 <= nx and nx < n) or c[ny][nx] == INT_MAX) {
                    int k = street_count;
                    street_count += 1;
                    dir ^= 1;  // flip the direction
                    ny = y;
                    nx = x;
                    bool found = false;
                    while (0 <= ny and ny < n and 0 <= nx and nx < n and c[ny][nx] != INT_MAX) {
                        if (intersection_index[ny][nx] != -1) {
                            streets_of[intersection_index[ny][nx]][dir == 1 ? 0 : 1] = k;
                            found = true;
                        }
                        ny += DIR_Y[dir];
                        nx += DIR_X[dir];
                    }
                    if (not found) {
                        street_count -= 1;
                    }
                }
            }
        }
    }
#ifdef VISUALIZE
    REP (i, intersection_count) {
        auto [y, x] = intersection_location[i];
        // cerr << "(" << y << ", " << x << ") looks " << streets_of[i][0] << " and " << streets_of[i][0] << endl;
    }
#endif
    REP (i, intersection_count) {
        assert (streets_of[i][0] != -1);
        assert (streets_of[i][1] != -1);
    }

    // Find start points of the intersections graph
    vector<int> start_intersections;
    if (intersection_index[start_y][start_x] != -1) {
        start_intersections.push_back(intersection_index[start_y][start_x]);
    } else {
        // use neighbors which is calculated before removing dead ends
        for (auto [ny, nx] : list_neighbors(start_y, start_x)) {
            int dy = ny - start_y;
            int dx = nx - start_x;
            while ((0 <= ny and ny < n and 0 <= nx and nx < n) and original_c[ny][nx] != INT_MAX) {
                if (intersection_index[ny][nx] != -1) {
                    start_intersections.push_back(intersection_index[ny][nx]);
                    break;
                }
                ny += dy;
                nx += dx;
            }
        }
    }

    auto make_random_path = [&]() -> pair<vector<int>, int> {
        vector<int> path;
        int score = 0;
        vector<bool> used(street_count);
        int cnt = 0;
        auto use = [&](int i) {
            score += (path.empty() ? dist_from_start_intersection(i) : dist[path.back()][i]);
            path.push_back(i);
            for (int k : streets_of[i]) {
                if (not used[k]) {
                    used[k] = true;
                    cnt += 1;
                }
            }
        };
        use(get_sample(start_intersections, gen));
        while (cnt < street_count) {
            int i = path.back();
            int j;
            while (true) {
                j = get_sample(g[i], gen).first;
                if (path.size() < 2 or path[path.size() - 2] != j or g[i].size() == 1) {
                    break;
                }
            }
            use(j);
        }
        score += dist_to_start_intersection(path.back());
        return {path, score};
    };

    vector<int> path;
    int score;
    tie(path, score) = make_random_path();

    vector<int> cumulative_score_left;
    vector<int> cumulative_score_right;
    vector<vector<bool>> cumulative_used_left;
    vector<vector<bool>> cumulative_used_right;
    vector<vector<int>> confluence_candidates(intersection_count);
    auto update_subinfo = [&]() {
        // left score
        cumulative_score_left.resize(path.size() + 1);
        REP (i, path.size()) {
            cumulative_score_left[i + 1] = cumulative_score_left[i] + (i == 0 ? dist_from_start_intersection(i) : dist[path[i - 1]][i]);
        }
        // right used
        cumulative_score_right.resize(path.size() + 1);
        REP_R (i, path.size()) {
            cumulative_score_right[i] = cumulative_score_right[i + 1] + (i == 0 ? dist_from_start_intersection(i) : dist[path[i - 1]][i]);
        }
        // left used
        cumulative_used_left.resize(path.size() + 1);
        cumulative_used_left[0] = vector<bool>(street_count);
        REP (i, path.size()) {
            cumulative_used_left[i + 1] = cumulative_used_left[i];
            for (int k : streets_of[path[i]]) {
                cumulative_used_left[i + 1][k] = true;
            }
        }
        // right used
        cumulative_used_right.resize(path.size() + 1);
        cumulative_used_right[path.size()] = vector<bool>(street_count);
        REP_R (i, path.size()) {
            cumulative_used_right[i] = cumulative_used_right[i + 1];
            for (int k : streets_of[path[i]]) {
                cumulative_used_right[i][k] = true;
            }
        }
        // confluence
        REP (k, street_count) {
            confluence_candidates[k].clear();
        }
        REP_R (i, n) {
            confluence_candidates[path[i]].push_back(i);
        }
    };
    update_subinfo();

    auto split_from = [&](const vector<int> &base_path, int split_start) -> tuple<vector<int>, int> {
        vector<int> path;
        int score = 0;
        vector<bool> used(street_count);
        int cnt = count(ALL(used), true);
        auto use = [&](int i) {
            score += (path.empty() ? dist_from_start_intersection(i) : dist[path.back()][i]);
            path.push_back(i);
            for (int k : streets_of[i]) {
                if (not used[k]) {
                    used[k] = true;
                    cnt += 1;
                }
            }
        };

        REP (i, split_start) {
            use(base_path[i]);
        }

        if (split_start == 0) {
            use(get_sample(start_intersections, gen));
        }
        while (cnt < street_count) {
            int i = (path.empty() ? base_path[split_start - 1] : path.back());
            int j;
            while (true) {
                j = get_sample(g[i], gen).first;
                if (path.size() < 2 or path[path.size() - 2] != j or g[i].size() == 1) {
                    break;
                }
            }
            use(j);

            // confluence
            for (int split_end : confluence_candidates[j]) {
                bool failed = true;
                REP (i, street_count) {
                    if (not (used[i] or cumulative_used_right[split_end][i])) {
                        failed = true;
                        break;
                    }
                }
                if (not failed) {
                    REP3 (k, split_end, base_path.size()) {
                        use(base_path[k]);
                        if (cnt == street_count) {
                            break;
                        }
                    }
                    assert (cnt == street_count);
                }
            }
        }
        score += dist_to_start_intersection(path.back());
        return {path, score};
    };

    vector<int> result = path;
    int highscore = score;

    int64_t iteration = 0;
    double temperature = 1.0;
    for (; ; ++ iteration) {
        if (iteration % 64 == 0) {
            chrono::high_resolution_clock::time_point clock_now = chrono::high_resolution_clock::now();
            temperature = static_cast<long double>((clock_end - clock_now).count()) / (clock_end - clock_begin).count();
            if (temperature <= 0.0) {
                cerr << "done  (iteration = " << iteration << ")" << endl;
                break;
            }
        }


        int split_start = uniform_int_distribution<int>(0, path.size() - 1)(gen);
        auto [next_path, next_score] = split_from(path, split_start);

        int delta = score - next_score;
        auto probability = [&]() {
            constexpr long double boltzmann = 0.05;
            return exp(boltzmann * delta / temperature);
        };
        if (delta >= 0 or bernoulli_distribution(probability())(gen)) {
            cerr << "update " << score << " -> " << next_score << endl;
            score = next_score;
            path = move(next_path);
            update_subinfo();

            if (score < highscore) {
                highscore = score;
                result = path;
            }
        }
    }

    // Reconstruct the path on the grid graph
    vector<pair<int, int>> reconstructed;
    reconstructed.emplace_back(start_y, start_x);
    for (int i : result) {
        reconstructed.push_back(intersection_location[i]);
    }
    while (true) {
        auto [y, x] = reconstructed.back();
        if (y == start_y and x == start_x) {
            break;
        }
        int i = intersection_index[y][x];
        assert (i != -1);
        if (find(ALL(start_intersections), i) != start_intersections.end()) {
            break;
        }
        int found = -1;
        for (auto [j, cost] : g[i]) {
            if (dist_to_start_intersection(j) + cost == dist_to_start_intersection(i)) {
                found = j;
                break;
            }
        }
        assert (found != -1);
        reconstructed.push_back(intersection_location[found]);
    }
    reconstructed.emplace_back(start_y, start_x);
#ifdef VISUALIZE
    for (auto [y, x] : reconstructed) {
        cerr << "(" << y << ", " << x << ") ";
    }
    cerr << endl;
#endif

    // Return the answer
    string ans = convert_to_command_string(reconstructed);
    // cerr << "ans = " << ans << endl;
    cerr << "score = " << highscore << endl;
    return ans;
}

int main() {
    constexpr auto TIME_LIMIT = chrono::milliseconds(3000);
    chrono::high_resolution_clock::time_point clock_begin = chrono::high_resolution_clock::now();
    xor_shift_128 gen(20210425);

    int n, sy, sx; cin >> n >> sy >> sx;
    vector<vector<int> > c(n, vector<int>(n));
    REP (y, n) {
        REP (x, n) {
            char a; cin >> a;
            assert (('5' <= a and a <= '9') or a == '#');
            c[y][x] = (a == '#' ? INT_MAX : a - '0');
        }
    }
    string ans = solve(n, sy, sx, c, gen, clock_begin + chrono::duration_cast<chrono::milliseconds>(TIME_LIMIT * 0.95));
    cout << ans << endl;
    return 0;
}
