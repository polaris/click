#include <alsa/asoundlib.h>

#include <boost/program_options.hpp>

#include <string>
#include <iostream>
#include <cmath>
#include <thread>
#include <atomic>
#include <chrono>

#include <signal.h>

bool verbose = false;
std::atomic<bool> running;

void initialize_alsa(snd_pcm_t** pcm, snd_pcm_status_t** status, snd_output_t** output, const std::string& device_name, unsigned int sample_rate, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels);
void set_hw_params(snd_pcm_t* pcm, unsigned int sample_rate, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels);
void set_sw_params(snd_pcm_t* pcm, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size);
void terminate_alsa(snd_pcm_t** pcm);
int playback_loop(snd_pcm_t* pcm, snd_pcm_status_t* status, snd_pcm_uframes_t period_size, unsigned int period_time, unsigned int sample_rate, unsigned int channels);
void generate_sine(const snd_pcm_channel_area_t *areas, snd_pcm_uframes_t offset, int count, double *_phase, unsigned int sample_rate, unsigned int channels);
int xrun_recovery(snd_pcm_t* pcm, int err);
void signal_handler(int);

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    std::string device_name = "default";
    unsigned int sample_rate = 48000;
    unsigned int period_time = 1000;
    unsigned int channels = 2;

    po::options_description desc("Options");
    desc.add_options()
        ("device,d", po::value<std::string>(&device_name)->default_value(device_name), "the device name of the audio hardware")
        ("samplerate,s", po::value<unsigned int>(&sample_rate)->default_value(sample_rate), "the sample rate in Hertz")
        ("periodtime,t", po::value<unsigned int>(&period_time)->default_value(period_time), "the packet time in us (125, 250, 333, 1000)")
        ("channels,c", po::value<unsigned int>(&channels)->default_value(channels), "the number of channels")
        ("verbose,v", "verbose output")
        ("help,h", "produce help message");

    try {
        po::variables_map vm;

        po::store(po::parse_command_line(argc, argv, desc), vm);

        po::notify(vm);
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }

        verbose = vm.count("verbose") > 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        std::cout << desc << "\n";
        return -1;
    }

    const auto period_size = static_cast<snd_pcm_uframes_t>(std::round(sample_rate * 0.000001 * period_time));
    const auto buffer_size = period_size * 2;

    snd_pcm_t* pcm = nullptr;
    snd_pcm_status_t* status = nullptr;
    snd_output_t* output = nullptr;

    initialize_alsa(&pcm, &status, &output, device_name, sample_rate, period_size, buffer_size, channels);

    if (verbose) {
        snd_pcm_dump(pcm, output);
    }

    running = true;

    signal(SIGINT, signal_handler);

    std::thread thread([pcm, status, period_size, period_time, sample_rate, channels] () {
        playback_loop(pcm, status, period_size, period_time, sample_rate, channels);
    });

    thread.join();

    terminate_alsa(&pcm);
}

void initialize_alsa(snd_pcm_t** pcm, snd_pcm_status_t** status, snd_output_t** output, const std::string& device_name, unsigned int sample_rate, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels) {
    int err = 0;
    
    err = snd_output_stdio_attach(output, stdout, 0);
    if (err < 0) {
        std::cerr << "Failed to attach stdio: " << snd_strerror(err) << "\n";
    }

    err = snd_pcm_open(pcm, device_name.c_str(), SND_PCM_STREAM_PLAYBACK, 0);
    if (err < 0) {
        std::cerr << "Failed to open device: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    snd_config_update_free_global();

    set_hw_params(*pcm, sample_rate, period_size, buffer_size, channels);
    
    set_sw_params(*pcm, period_size, buffer_size);

    snd_pcm_status_alloca(status);
}

void set_hw_params(snd_pcm_t* pcm, unsigned int sample_rate, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels) {
    int err = 0;

    snd_pcm_hw_params_t* params = nullptr;
    err = snd_pcm_hw_params_malloc(&params);
    if (err < 0) {
        std::cerr << "Failed to allocate snd_pcm_hw_params_t: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_any(pcm, params);
    if (err < 0) {
        std::cerr << "Failed to fill hardware params with configuration space for a PCM: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_access(pcm, params, SND_PCM_ACCESS_MMAP_INTERLEAVED);
    if (err < 0) {
        std::cerr << "Failed to set access type to SND_PCM_ACCESS_MMAP_INTERLEAVED: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_format(pcm, params, SND_PCM_FORMAT_S16_LE);
    if (err < 0) {
        std::cerr << "Failed to set format to SND_PCM_FORMAT_S16_LE: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_rate(pcm, params, sample_rate, 0);
    if (err < 0) {
        std::cerr << "Failed to set sample rate to " << sample_rate << ": " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_channels(pcm, params, channels);
    if (err < 0) {
        std::cerr << "Failed to set number of channels to " << channels << " : " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_period_size(pcm, params, period_size, 0);
    if (err < 0) {
        std::cerr << "Failed to set period size to " << period_size << " : " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_periods(pcm, params, 2, 0);
    if (err < 0) {
        std::cerr << "Failed to set periods to 2: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params_set_buffer_size(pcm, params, buffer_size);
    if (err < 0) {
        std::cerr << "Failed to set buffer size to " << buffer_size << ": " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_hw_params(pcm, params);
    if (err < 0) {
        std::cerr << "Failed to install hardware configuration for PCM: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    snd_pcm_hw_params_free(params);
}

void set_sw_params(snd_pcm_t* pcm, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size) {
    int err = 0;

    snd_pcm_sw_params_t* params = nullptr;
    err = snd_pcm_sw_params_malloc(&params);
    if (err < 0) {
        std::cerr << "Failed to allocate snd_pcm_sw_params_t: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_sw_params_current(pcm, params);
    if (err < 0) {
        std::cerr << "Failed to get current software configuration: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_sw_params_set_avail_min(pcm, params, period_size);
    if (err < 0) {
        std::cerr << "Failed to set avail min: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_sw_params_set_start_threshold(pcm, params, (buffer_size / period_size) * period_size);
    if (err < 0) {
        std::cerr << "Failed to set start threshold: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_sw_params_set_tstamp_mode(pcm, params, SND_PCM_TSTAMP_ENABLE);
    if (err < 0) {
        std::cerr << "Failed to set tstamp mode to SND_PCM_TSTAMP_ENABLE: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    err = snd_pcm_sw_params(pcm, params);
    if (err < 0) {
        std::cerr << "Failed to install software configuration for PCM: " << snd_strerror(err) << "\n";
        exit(-1);
    }

    snd_pcm_sw_params_free(params);
}

int playback_loop(snd_pcm_t* pcm, snd_pcm_status_t* status, snd_pcm_uframes_t period_size, unsigned int period_time, unsigned int sample_rate, unsigned int channels) {
    int err = 0, first = 1;
    double phase = 0;

    while (running) {
        auto state = snd_pcm_state(pcm);
        if (state == SND_PCM_STATE_XRUN) {
            err = xrun_recovery(pcm, -EPIPE);
            if (err < 0) {
                std::cerr << "Failed to recover from XRUN: " << snd_strerror(err) << "\n";
                return err;
            }
            first = 1;
        } else if (state == SND_PCM_STATE_SUSPENDED) {
            err = xrun_recovery(pcm, -ESTRPIPE);
            if (err < 0) {
                std::cerr << "Failed to recover from SUSPEND: " << snd_strerror(err) << "\n";
                return err;
            }
        }
        snd_pcm_sframes_t avail = snd_pcm_avail_update(pcm);
        if (avail < 0) {
                err = xrun_recovery(pcm, static_cast<int>(avail));
                if (err < 0) {
                    std::cerr << "Failed to update avail: " << snd_strerror(err) << "\n";
                    return err;
                }
                first = 1;
                continue;
        }
        if (avail < static_cast<snd_pcm_sframes_t>(period_size)) {
            if (first) {
                first = 0;
                err = snd_pcm_start(pcm);
                if (err < 0) {
                    std::cerr << "Failed to start: " << snd_strerror(err) << "\n";
                    return err;
                }
            } else {
                err = snd_pcm_wait(pcm, -1);
                if (err < 0) {
                    if ((err = xrun_recovery(pcm, err)) < 0) {
                        std::cerr << "Failed to wait for PCM: " << snd_strerror(err) << "\n";
                        return err;
                    }
                    first = 1;
                }
            }
            continue;
        }

        err = snd_pcm_status(pcm, status);
        assert(err == 0);

        struct timespec timestamp;
        snd_pcm_status_get_htstamp(status, &timestamp);
        snd_pcm_sframes_t delay = snd_pcm_status_get_delay(status);
        state = snd_pcm_status_get_state(status);

        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        long us = static_cast<long>(std::round(static_cast<double>(now.tv_nsec) / 1000.0));
        long us2 = static_cast<long>(std::round(static_cast<double>(timestamp.tv_nsec) / 1000.0));
        long diff = 1000000 - us;
        long diff2 = 1000000 - us2;
        if (diff < period_time) {
            printf("%03ld, %03ld, %ld, %ld\n", diff, diff2, avail, delay);
        }

        snd_pcm_uframes_t size = period_size;
        while (size > 0) {
            snd_pcm_uframes_t offset = 0, frames = size;
            const snd_pcm_channel_area_t* channel_area = nullptr;
            err = snd_pcm_mmap_begin(pcm, &channel_area, &offset, &frames);
            if (err < 0) {
                if ((err = xrun_recovery(pcm, err)) < 0) {
                    std::cerr << "Failed MMAP begin: " << snd_strerror(err) << "\n";
                    return err;
                }
                first = 1;
            }
            generate_sine(channel_area, offset, static_cast<int>(frames), &phase, sample_rate, channels);
            snd_pcm_sframes_t commit_result = snd_pcm_mmap_commit(pcm, offset, frames);
            if (commit_result < 0 || (snd_pcm_uframes_t)commit_result != frames) {
                if ((err = xrun_recovery(pcm, commit_result >= 0 ? -EPIPE : static_cast<int>(commit_result))) < 0) {
                    std::cerr << "Failed MMAP commit: " << snd_strerror(err) << "\n";
                    return err;
                }
                first = 1;
            }
            size -= frames;
        }
    }
    return 0;
}

void generate_sine(const snd_pcm_channel_area_t *areas, snd_pcm_uframes_t offset, int count, double *_phase, unsigned int sample_rate, unsigned int channels) {
    static double freq = 440;                               /* sinusoidal wave frequency in Hz */
    static double max_phase = 2. * M_PI;
    double phase = *_phase;
    double step = max_phase*freq/(double)sample_rate;
    unsigned char** samples = new unsigned char* [channels];
    int* steps = new int [channels];
    unsigned int chn;
    int format_bits = snd_pcm_format_width(SND_PCM_FORMAT_S16_LE);
    unsigned int maxval = (1 << (format_bits - 1)) - 1;
    int bps = format_bits / 8;  /* bytes per sample */
    int phys_bps = snd_pcm_format_physical_width(SND_PCM_FORMAT_S16_LE) / 8;
    int big_endian = snd_pcm_format_big_endian(SND_PCM_FORMAT_S16_LE) == 1;
    int to_unsigned = snd_pcm_format_unsigned(SND_PCM_FORMAT_S16_LE) == 1;
    int is_float = (SND_PCM_FORMAT_S16_LE == SND_PCM_FORMAT_FLOAT_LE ||
                    SND_PCM_FORMAT_S16_LE == SND_PCM_FORMAT_FLOAT_BE);
    /* verify and prepare the contents of areas */
    for (chn = 0; chn < channels; chn++) {
            if ((areas[chn].first % 8) != 0) {
                    printf("areas[%i].first == %i, aborting...\n", chn, areas[chn].first);
                    exit(EXIT_FAILURE);
            }
            samples[chn] = /*(signed short *)*/(((unsigned char *)areas[chn].addr) + (areas[chn].first / 8));
            if ((areas[chn].step % 16) != 0) {
                    printf("areas[%i].step == %i, aborting...\n", chn, areas[chn].step);
                    exit(EXIT_FAILURE);
            }
            steps[chn] = areas[chn].step / 8;
            samples[chn] += offset * steps[chn];
    }
    /* fill the channel areas */
    while (count-- > 0) {
            union {
                    float f;
                    int i;
            } fval;
            int res, i;
            if (is_float) {
                    fval.f = sin(phase);
                    res = fval.i;
            } else
                    res = sin(phase) * maxval;
            if (to_unsigned)
                    res ^= 1U << (format_bits - 1);
            for (chn = 0; chn < channels; chn++) {
                    /* Generate data in native endian format */
                    if (big_endian) {
                            for (i = 0; i < bps; i++)
                                    *(samples[chn] + phys_bps - 1 - i) = (res >> i * 8) & 0xff;
                    } else {
                            for (i = 0; i < bps; i++)
                                    *(samples[chn] + i) = (res >>  i * 8) & 0xff;
                    }
                    samples[chn] += steps[chn];
            }
            phase += step;
            if (phase >= max_phase)
                    phase -= max_phase;
    }
    *_phase = phase;

    delete [] samples;
    delete [] steps;
}

void terminate_alsa(snd_pcm_t** pcm) {
    if (*pcm != nullptr) {
        const int err = snd_pcm_close(*pcm);
        if (err < 0) {
            std::cerr << "Failed to close PCM: " << snd_strerror(err) << "\n";
        }
        *pcm = nullptr;
    }
}

int xrun_recovery(snd_pcm_t* pcm, int err) {
    if (verbose) {
        std::cerr << "Stream recovery\n";
    }
    if (err == -EPIPE) {
        err = snd_pcm_prepare(pcm);
        if (err < 0) {
            std::cerr << "Can't recovery from underrun, prepare failed: " << snd_strerror(err) << "\n";
            return err;
        }
        return 0;
    } else if (err == -ESTRPIPE) {
        while ((err = snd_pcm_resume(pcm)) == -EAGAIN) {
            sleep(1); // Wait until the suspend flag is released
        }
        if (err < 0) {
            err = snd_pcm_prepare(pcm);
            if (err < 0) {
                std::cerr << "Can't recovery from suspend, prepare failed: " << snd_strerror(err) << "\n";
                return err;
            }
        }
        return 0;
    }
    return err;
}

void signal_handler(int) {
    if (verbose) {
        std::cerr << "Shutting down\n";
    }
    running = false;
}
