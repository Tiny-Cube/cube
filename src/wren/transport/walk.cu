// One CUDA thread follows a whole independent continuous-time trajectory.
// State words are interleaved by lane for coalesced reads across trajectories.
__device__ unsigned int mt_word(unsigned int* state, int width, int& cursor) {
    const int next = (cursor + 1) % 624;
    const unsigned int mixed = (state[cursor * width] & 0x80000000u)
                             | (state[next * width] & 0x7fffffffu);
    unsigned int word = state[((cursor + 397) % 624) * width]
                      ^ (mixed >> 1) ^ ((mixed & 1u) ? 0x9908b0dfu : 0u);
    state[cursor * width] = word;
    cursor = next;
    word ^= word >> 11;
    word ^= (word << 7) & 0x9d2c5680u;
    word ^= (word << 15) & 0xefc60000u;
    return word ^ (word >> 18);
}

__device__ double positive_uniform(unsigned int* state, int width, int& cursor) {
    double u;
    do {
        const unsigned int high = mt_word(state, width, cursor) >> 5;
        const unsigned int low = mt_word(state, width, cursor) >> 6;
        u = (high * 67108864.0 + low) / 9007199254740992.0;
    } while (u == 0.0);
    return u;
}

extern "C" __global__ void kmc_walk(
    int width, const int* row, const int* target, const int* event_id,
    const double* cdf, const double* delta, const int* initial, unsigned int* state,
    int nsamples, const double* grid, long long max_events, double* observations,
    long long* counts, int* capped, int* errors, int trace,
    const long long* trace_offset, double* times, int* events
) {
    const int lane = blockIdx.x * blockDim.x + threadIdx.x;
    if (lane >= width) return;
    unsigned int* rng = state + lane;
    int cursor = 0;
    int site = initial[lane];
    int sample = 0;
    double t = 0.0;
    double xyz[3] = {0.0, 0.0, 0.0};
    long long steps = 0;
    const double end = grid[nsamples - 1];
    const long long sample_start = (long long)lane * nsamples * 3;
    long long entry = trace ? trace_offset[lane] : 0;
    while (t < end && steps < max_events) {
        const int begin = row[site], stop = row[site + 1];
        int edge = -1;
        double next = end;
        if (begin < stop) {
            const double total = cdf[stop - 1];
            const double dt = -log(positive_uniform(rng, width, cursor)) / total;
            if (dt <= end - t) {
                next = t + dt;
                if (!isfinite(next) || next <= t) { errors[lane] = 1; return; }
                const double threshold = positive_uniform(rng, width, cursor) * total;
                // First cumulative rate strictly above threshold; protect the
                // last positive channel when threshold rounds up to total.
                int lo = begin, hi = stop - 1;
                while (lo < hi) {
                    int mid = lo + (hi - lo) / 2;
                    if (cdf[mid] <= threshold) lo = mid + 1; else hi = mid;
                }
                edge = lo;
            }
        }
        if (!trace) {
            while (sample < nsamples && grid[sample] < next) {
                for (int a = 0; a < 3; ++a) observations[sample_start + 3*sample+a] = xyz[a];
                ++sample;
            }
        } else {
            if (entry >= trace_offset[lane+1]) { errors[lane] = 2; return; }
            times[entry] = next;
            events[entry++] = edge >= 0 ? event_id[edge] : -1;
        }
        t = next;
        if (edge < 0) break;
        for (int a = 0; a < 3; ++a) xyz[a] += delta[(long long)edge * 3 + a];
        site = target[edge];
        ++steps;
    }
    const int limited = t < end;
    counts[lane] = steps;
    capped[lane] = limited;
    if (!trace) {
        for (; sample < nsamples; ++sample) {
            for (int a = 0; a < 3; ++a)
                observations[sample_start + 3*sample+a] = limited ? nan("") : xyz[a];
        }
    }
}
