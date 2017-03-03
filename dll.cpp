#include <alsa/asoundlib.h>

#include <boost/program_options.hpp>

#include <string>
#include <iostream>
#include <cmath>
#include <thread>
#include <atomic>
#include <chrono>

#include <signal.h>
#include <inttypes.h>

bool verbose = false;
std::atomic<bool> running;

void initialize_alsa(snd_pcm_t** pcm, snd_output_t** output, const std::string& device_name, unsigned int sample_rate,
    snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels);
void set_hw_params(snd_pcm_t* pcm, unsigned int sample_rate, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels);
void set_sw_params(snd_pcm_t* pcm, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size);
void terminate_alsa(snd_pcm_t** pcm);
int playback_loop(snd_pcm_t* pcm, snd_pcm_uframes_t period_size, unsigned int period_time_us, unsigned int sample_rate, unsigned int channels);
void generate_click_track(int64_t system_time_us, unsigned int click_duration_us, unsigned int sample_rate,
    unsigned int channels, double* phase, const snd_pcm_channel_area_t* channel_area, snd_pcm_uframes_t offset, snd_pcm_uframes_t frames);
int xrun_recovery(snd_pcm_t* pcm, int err);
void signal_handler(int);
uint64_t timespec_us(const struct timespec *ts);

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    std::string device_name = "default";
    unsigned int sample_rate = 48000;
    unsigned int period_time_us = 1000;
    unsigned int channels = 2;

    po::options_description desc("Options");
    desc.add_options()
        ("device,d", po::value<std::string>(&device_name)->default_value(device_name), "the device name of the audio hardware")
        ("samplerate,s", po::value<unsigned int>(&sample_rate)->default_value(sample_rate), "the sample rate in Hertz")
        ("periodtime,t", po::value<unsigned int>(&period_time_us)->default_value(period_time_us), "the packet time in us (125, 250, 333, 1000)")
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

    const auto period_size = static_cast<snd_pcm_uframes_t>(std::round(sample_rate * 0.000001 * period_time_us));
    const auto buffer_size = period_size * 2;

    snd_pcm_t* pcm = nullptr;
    snd_output_t* output = nullptr;

    initialize_alsa(&pcm, &output, device_name, sample_rate, period_size, buffer_size, channels);

    if (verbose) {
        snd_pcm_dump(pcm, output);
    }

    running = true;

    signal(SIGINT, signal_handler);

    std::thread thread([pcm, period_size, period_time_us, sample_rate, channels] () {
        playback_loop(pcm, period_size, period_time_us, sample_rate, channels);
    });

    thread.join();

    terminate_alsa(&pcm);
}

void initialize_alsa(snd_pcm_t** pcm, snd_output_t** output, const std::string& device_name, unsigned int sample_rate, snd_pcm_uframes_t period_size, snd_pcm_uframes_t buffer_size, unsigned int channels) {
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

    int supports_audio_wallclock_ts = snd_pcm_hw_params_supports_audio_wallclock_ts(params);
    if (supports_audio_wallclock_ts == 1) {
        std::cout << "Audio wallclock timestamps supported\n";
    } else {
        std::cout << "Audio wallclock timestamps not supported\n";
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

int64_t frames_us(unsigned int sample_rate, snd_pcm_sframes_t frames) {
    return static_cast<int64_t>((1.0 / sample_rate) * static_cast<double>(1000000 * frames));
}

double read_timer() {
    struct timespec realtime;
    clock_gettime(CLOCK_REALTIME, &realtime);
    return static_cast<double>(timespec_us(&realtime));
}

int playback_loop(snd_pcm_t* pcm, snd_pcm_uframes_t period_size, unsigned int period_time_us, unsigned int sample_rate, unsigned int channels) {
    int err = 0, first = 1;
    double phase = 0;

    double sqrt2 = 1.414213562373095;
    double pi    = 3.141592653589793;

    double tper = period_time_us;
    double omega = 2 * pi * 1.25e-7 * period_time_us;
    double a = 0;
    double b = sqrt2 * omega;
    double c = omega * omega;
    double e2 = tper;
    double t0 = read_timer();
    double t1 = t0 + e2;
    unsigned int cnt = 0;

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

        snd_pcm_sframes_t avail = 0, delay = 0;
        err = snd_pcm_avail_delay(pcm, &avail, &delay);
        if (err < 0) {
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

        struct timespec system_time;
        clock_gettime(CLOCK_REALTIME, &system_time);
        const int64_t system_time_us = timespec_us(&system_time);

        double e = static_cast<double>(system_time_us) - t1;
        t0 = t1;
        t1 += b * e + e2;
        e2 += c * e;
        cnt += 1;
        if (cnt % 1000 == 0) {
            std::cout.precision(32);
            std::cout << t1 - t0 << "\n"; // << ", " << t1 << std::endl;
        }

        const int64_t delay_us = frames_us(sample_rate, delay);
        int64_t presentation_time_us = system_time_us + delay_us;

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

            generate_click_track(presentation_time_us, 20000, sample_rate, channels, &phase, channel_area, offset, frames);

            presentation_time_us += frames_us(sample_rate, frames);

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

void generate_click_track(int64_t presentation_time_us, unsigned int click_duration_us, 
                            unsigned int sample_rate, unsigned int channels, double* phase, 
                            const snd_pcm_channel_area_t* channel_area, snd_pcm_uframes_t offset, snd_pcm_uframes_t frames) {
    static const double freq = 880;
    static const double max_phase = 2. * M_PI;

    double phs = *phase;

    const double step = max_phase * freq / static_cast<double>(sample_rate);
    const double inc = 1000.0 / sample_rate;

    const int64_t floor_second = 1000000 * (presentation_time_us / 1000000);
    double diff_second = static_cast<double>(presentation_time_us - floor_second);

    for (unsigned int frame = 0; frame < frames; frame++) {
        int16_t value = 0;
        if (diff_second < click_duration_us) {
            value = static_cast<int16_t>(sin(phs) * (.75 * 0x8000));
            phs += step;
            if (phs >= max_phase) {
                phs -= max_phase;
            }
        } else {
            phs = 0;
        }
        diff_second += inc;

        for (unsigned int channel = 0; channel < channels; channel++) {
            auto sampleAddress = (int16_t*)channel_area[channel].addr;
            sampleAddress += (channel_area[channel].first + (frame + offset) * channel_area[channel].step) / 16;
            *sampleAddress = value;
        }
    }
    *phase = phs;
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

uint64_t timespec_us(const struct timespec *ts) {
    return ts->tv_sec * 1000000LLU + ts->tv_nsec / 1000LLU;
}
