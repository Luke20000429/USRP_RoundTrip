//
// Copyright 2014 Ettus Research LLC
// Copyright 2018 Ettus Research, a National Instruments Company
//
// SPDX-License-Identifier: GPL-3.0-or-later
//

// This tiny program is meant as an example on how to set up UHD
// projects using CMake.
// The program itself only initializes a USRP. For more elaborate examples,
// have a look at the files in host/examples/.

#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/static.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <csignal>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>
#include <chrono>

namespace po = boost::program_options;

#define NOW() (std::chrono::high_resolution_clock::now())
using timestamp_t = std::chrono::high_resolution_clock::time_point;

float zc_real[] = {1.        , 0.99054969, 0.91601557, 0.67816903, 0.19368445,-0.47329161,
 -0.96835464,-0.75779363, 0.23845923, 0.99579618, 0.28273252,-0.9414233 ,
 -0.2606644 , 0.99894849,-0.30465196,-0.69484163, 0.99054969,-0.58963889,
 -0.03439022, 0.53272534,-0.81437437, 0.93344339,-0.96835464, 0.96237739,
 -0.9065768 , 0.74263314,-0.39061386,-0.1711384 , 0.77255565,-0.98714441,
  0.41161885, 0.643763  ,-0.92497267,-0.21612866, 0.98322007, 0.28273252,
 -0.86410814,-0.78691144, 0.14850236, 0.87542143, 0.96237739, 0.60800312,
  0.14850236,-0.21612866,-0.43240741,-0.51318047,-0.47329161,-0.30465196,
  0.01146542, 0.4529686 , 0.87542143, 0.97382271, 0.4529686 ,-0.51318047,
 -0.99343412,-0.1711384 , 0.94890818, 0.32641121,-0.97877872, 0.05729694,
  0.89666134,-0.84012465, 0.10300796, 0.60800312,-0.9558941 , 0.98322007,
 -0.86410814, 0.74263314,-0.69484163, 0.74263314,-0.86410814, 0.98322007,
 -0.9558941 , 0.60800312, 0.10300796,-0.84012465, 0.89666134, 0.05729694,
 -0.97877872, 0.32641121, 0.94890818,-0.1711384 ,-0.99343412,-0.51318047,
  0.4529686 , 0.97382271, 0.87542143, 0.4529686 , 0.01146542,-0.30465196,
 -0.47329161,-0.51318047,-0.43240741,-0.21612866, 0.14850236, 0.60800312,
  0.96237739, 0.87542143, 0.14850236,-0.78691144,-0.86410814, 0.28273252,
  0.98322007,-0.21612866,-0.92497267, 0.643763  , 0.41161885,-0.98714441,
  0.77255565,-0.1711384 ,-0.39061386, 0.74263314,-0.9065768 , 0.96237739,
 -0.96835464, 0.93344339,-0.81437437, 0.53272534,-0.03439022,-0.58963889,
  0.99054969,-0.69484163,-0.30465196, 0.99894849,-0.2606644 ,-0.9414233 ,
  0.28273252, 0.99579618, 0.23845923,-0.75779363,-0.96835464,-0.47329161,
  0.19368445, 0.67816903, 0.91601557, 0.99054969, 1.        };

float zc_imag[] = {0.00000000e+00,-1.37154309e-01,-4.01142723e-01,-7.34905956e-01,
 -9.81063878e-01,-8.80905813e-01,-2.49578219e-01, 6.52494307e-01,
  9.71152510e-01, 9.15967646e-02,-9.59198793e-01,-3.37227189e-01,
  9.65429474e-01, 4.58465904e-02,-9.52463743e-01, 7.19162783e-01,
  1.37154309e-01,-8.07666999e-01, 9.99408482e-01,-8.46288194e-01,
  5.80339894e-01,-3.58724727e-01, 2.49578219e-01,-2.71716318e-01,
  4.22040875e-01,-6.69698454e-01, 9.20554624e-01,-9.85246999e-01,
  6.34947064e-01, 1.59830884e-01,-9.11356088e-01, 7.65224934e-01,
  3.80033640e-01,-9.76364892e-01,-1.82423417e-01, 9.59198793e-01,
  5.03306193e-01,-6.17065950e-01,-9.88912054e-01,-4.83360453e-01,
  2.71716318e-01, 7.93934633e-01, 9.88912054e-01, 9.76364892e-01,
  9.01678340e-01, 8.58280727e-01, 8.80905813e-01, 9.52463743e-01,
  9.99934270e-01, 8.91526469e-01, 4.83360453e-01,-2.27308885e-01,
 -8.91526469e-01,-8.58280727e-01, 1.14405616e-01, 9.85246999e-01,
  3.15552329e-01,-9.45227866e-01,-2.04920027e-01, 9.98357181e-01,
 -4.42717108e-01,-5.42393375e-01, 9.94680532e-01,-7.93934633e-01,
  2.93711544e-01, 1.82423417e-01,-5.03306193e-01, 6.69698454e-01,
 -7.19162783e-01, 6.69698454e-01,-5.03306193e-01, 1.82423417e-01,
  2.93711544e-01,-7.93934633e-01, 9.94680532e-01,-5.42393375e-01,
 -4.42717108e-01, 9.98357181e-01,-2.04920027e-01,-9.45227866e-01,
  3.15552329e-01, 9.85246999e-01, 1.14405616e-01,-8.58280727e-01,
 -8.91526469e-01,-2.27308885e-01, 4.83360453e-01, 8.91526469e-01,
  9.99934270e-01, 9.52463743e-01, 8.80905813e-01, 8.58280727e-01,
  9.01678340e-01, 9.76364892e-01, 9.88912054e-01, 7.93934633e-01,
  2.71716318e-01,-4.83360453e-01,-9.88912054e-01,-6.17065950e-01,
  5.03306193e-01, 9.59198793e-01,-1.82423417e-01,-9.76364892e-01,
  3.80033640e-01, 7.65224934e-01,-9.11356088e-01, 1.59830884e-01,
  6.34947064e-01,-9.85246999e-01, 9.20554624e-01,-6.69698454e-01,
  4.22040875e-01,-2.71716318e-01, 2.49578219e-01,-3.58724727e-01,
  5.80339894e-01,-8.46288194e-01, 9.99408482e-01,-8.07666999e-01,
  1.37154309e-01, 7.19162783e-01,-9.52463743e-01, 4.58465904e-02,
  9.65429474e-01,-3.37227189e-01,-9.59198793e-01, 9.15967646e-02,
  9.71152510e-01, 6.52494307e-01,-2.49578219e-01,-8.80905813e-01,
 -9.81063878e-01,-7.34905956e-01,-4.01142723e-01,-1.37154309e-01,
  7.83872988e-14};

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int)
{
    stop_signal_called = true;
}


int UHD_SAFE_MAIN(int argc, char* argv[])
{
    // transmit variables to be set by po
    std::string tx_args, wave_type, tx_ant, tx_subdev, ref, otw, tx_channels;
    double tx_rate, tx_freq, tx_gain, wave_freq, tx_bw;
    // float ampl;

    // receive variables to be set by po
    std::string rx_args, file, type, rx_ant, rx_subdev, rx_channels;
    size_t total_num_samps, spb;
    double rx_rate, rx_freq, rx_gain, rx_bw;
    double settling;

    // setup the program options
    po::options_description desc("Allowed options");
    // clang-format off
    desc.add_options()
        ("help", "help message")
        ("tx-args", po::value<std::string>(&tx_args)->default_value(""), "uhd transmit device address args")
        ("rx-args", po::value<std::string>(&rx_args)->default_value(""), "uhd receive device address args")
        // ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "sample type in file: double, float, or short")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("settling", po::value<double>(&settling)->default_value(double(0.2)), "settling time (seconds) before receiving")
        ("spb", po::value<size_t>(&spb)->default_value(8), "samples per buffer, 0 for default")
        ("tx-rate", po::value<double>(&tx_rate)->default_value(double(40.0e6)), "rate of transmit outgoing samples, 20MHz by default")
        ("rx-rate", po::value<double>(&rx_rate)->default_value(double(40.0e6)), "rate of receive incoming samples, 20MHz by default")
        ("tx-freq", po::value<double>(&tx_freq)->default_value(double(2.45e9)), "transmit RF center frequency in Hz, 2.45GHz by default")
        ("rx-freq", po::value<double>(&rx_freq)->default_value(double(5.0e9)), "receive RF center frequency in Hz, 5GHz by default")
        // ("ampl", po::value<float>(&ampl)->default_value(float(0.3)), "amplitude of the waveform [0 to 0.7]")
        ("tx-gain", po::value<double>(&tx_gain)->default_value(double(90.0)), "gain for the transmit RF chain, 70dB by default")
        ("rx-gain", po::value<double>(&rx_gain)->default_value(double(40.0)), "gain for the receive RF chain, 40dB by default")
        ("tx-ant", po::value<std::string>(&tx_ant), "transmit antenna selection")
        ("rx-ant", po::value<std::string>(&rx_ant), "receive antenna selection")
        ("tx-subdev", po::value<std::string>(&tx_subdev), "transmit subdevice specification")
        ("rx-subdev", po::value<std::string>(&rx_subdev), "receive subdevice specification")
        ("tx-bw", po::value<double>(&tx_bw)->default_value(double(20.0e6)), "analog transmit filter bandwidth in Hz")
        ("rx-bw", po::value<double>(&rx_bw)->default_value(double(20.0e6)), "analog receive filter bandwidth in Hz")
        // ("wave-type", po::value<std::string>(&wave_type)->default_value("SINE"), "waveform type (CONST, SQUARE, RAMP, SINE)")
        // ("wave-freq", po::value<double>(&wave_freq)->default_value(1000), "waveform frequency in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "clock reference (internal, external, mimo)")
        ("otw", po::value<std::string>(&otw)->default_value("sc16"), "specify the over-the-wire sample mode")
        ("tx-channels", po::value<std::string>(&tx_channels)->default_value("0"), "which TX channel(s) to use (specify \"0\", \"1\", \"0,1\", etc), 0 by default")
        ("rx-channels", po::value<std::string>(&rx_channels)->default_value("1"), "which RX channel(s) to use (specify \"0\", \"1\", \"0,1\", etc), 1 by default")
        ("tx-int-n", "tune USRP TX with integer-N tuning")
        ("rx-int-n", "tune USRP RX with integer-N tuning")
    ;
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // print the help message
    if (vm.count("help")) {
        std::cout << "UHD TXRX Loopback to File " << desc << std::endl;
        return ~0;
    }

    // create a usrp device
    std::cout << boost::format("Creating the transmit usrp device with: %s...") % tx_args
              << std::endl;
    uhd::usrp::multi_usrp::sptr tx_usrp = uhd::usrp::multi_usrp::make(tx_args);
    std::cout << std::endl;
    std::cout << boost::format("Creating the receive usrp device with: %s...") % rx_args
              << std::endl;
    uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

    // always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("tx-subdev"))
        tx_usrp->set_tx_subdev_spec(tx_subdev);
    if (vm.count("rx-subdev"))
        rx_usrp->set_rx_subdev_spec(rx_subdev);

    // detect which channels to use
    std::vector<std::string> tx_channel_strings;
    std::vector<size_t> tx_channel_nums;
    boost::split(tx_channel_strings, tx_channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < tx_channel_strings.size(); ch++) {
        size_t chan = std::stoi(tx_channel_strings[ch]);
        if (chan >= tx_usrp->get_tx_num_channels()) {
            throw std::runtime_error("Invalid TX channel(s) specified.");
        } else
            tx_channel_nums.push_back(std::stoi(tx_channel_strings[ch]));
    }
    std::vector<std::string> rx_channel_strings;
    std::vector<size_t> rx_channel_nums;
    boost::split(rx_channel_strings, rx_channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < rx_channel_strings.size(); ch++) {
        size_t chan = std::stoi(rx_channel_strings[ch]);
        if (chan >= rx_usrp->get_rx_num_channels()) {
            throw std::runtime_error("Invalid RX channel(s) specified.");
        } else
            rx_channel_nums.push_back(std::stoi(rx_channel_strings[ch]));
    }

    // Lock mboard clocks
    if (vm.count("ref")) {
        tx_usrp->set_clock_source(ref);
        rx_usrp->set_clock_source(ref);
    }

    std::cout << "Using TX Device: " << tx_usrp->get_pp_string() << std::endl;
    std::cout << "Using RX Device: " << rx_usrp->get_pp_string() << std::endl;

    // set the transmit sample rate
    if (not vm.count("tx-rate")) {
        std::cerr << "Please specify the transmit sample rate with --tx-rate"
                  << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (tx_rate / 1e6)
              << std::endl;
    tx_usrp->set_tx_rate(tx_rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...")
                     % (tx_usrp->get_tx_rate() / 1e6)
              << std::endl
              << std::endl;

    // set the receive sample rate
    if (not vm.count("rx-rate")) {
        std::cerr << "Please specify the sample rate with --rx-rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rx_rate / 1e6)
              << std::endl;
    rx_usrp->set_rx_rate(rx_rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...")
                     % (rx_usrp->get_rx_rate() / 1e6)
              << std::endl
              << std::endl;

    // set the transmit center frequency
    if (not vm.count("tx-freq")) {
        std::cerr << "Please specify the transmit center frequency with --tx-freq"
                  << std::endl;
        return ~0;
    }

    for (size_t ch = 0; ch < tx_channel_nums.size(); ch++) {
        size_t channel = tx_channel_nums[ch];
        if (tx_channel_nums.size() > 1) {
            std::cout << "Configuring TX Channel " << channel << std::endl;
        }
        std::cout << boost::format("Setting TX Freq: %f MHz...") % (tx_freq / 1e6)
                  << std::endl;
        uhd::tune_request_t tx_tune_request(tx_freq);
        if (vm.count("tx-int-n"))
            tx_tune_request.args = uhd::device_addr_t("mode_n=integer");
        tx_usrp->set_tx_freq(tx_tune_request, channel);
        std::cout << boost::format("Actual TX Freq: %f MHz...")
                         % (tx_usrp->get_tx_freq(channel) / 1e6)
                  << std::endl
                  << std::endl;

        // set the rf gain
        if (vm.count("tx-gain")) {
            std::cout << boost::format("Setting TX Gain: %f dB...") % tx_gain
                      << std::endl;
            // tx_usrp->set_tx_gain(tx_gain, channel);
            tx_usrp->set_normalized_tx_gain(1.0, channel);
            std::cout << boost::format("Actual TX Gain: %f dB...")
                             % tx_usrp->get_tx_gain(channel)
                      << std::endl
                      << std::endl;
        }

        // set the analog frontend filter bandwidth
        if (vm.count("tx-bw")) {
            std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % tx_bw
                      << std::endl;
            tx_usrp->set_tx_bandwidth(tx_bw, channel);
            std::cout << boost::format("Actual TX Bandwidth: %f MHz...")
                             % tx_usrp->get_tx_bandwidth(channel)
                      << std::endl
                      << std::endl;
        }

        // set the antenna
        if (vm.count("tx-ant"))
            tx_usrp->set_tx_antenna(tx_ant, channel);
    }

    for (size_t ch = 0; ch < rx_channel_nums.size(); ch++) {
        size_t channel = rx_channel_nums[ch];
        if (rx_channel_nums.size() > 1) {
            std::cout << "Configuring RX Channel " << channel << std::endl;
        }

        // set the receive center frequency
        if (not vm.count("rx-freq")) {
            std::cerr << "Please specify the center frequency with --rx-freq"
                      << std::endl;
            return ~0;
        }
        std::cout << boost::format("Setting RX Freq: %f MHz...") % (rx_freq / 1e6)
                  << std::endl;
        uhd::tune_request_t rx_tune_request(rx_freq);
        if (vm.count("rx-int-n"))
            rx_tune_request.args = uhd::device_addr_t("mode_n=integer");
        rx_usrp->set_rx_freq(rx_tune_request, channel);
        std::cout << boost::format("Actual RX Freq: %f MHz...")
                         % (rx_usrp->get_rx_freq(channel) / 1e6)
                  << std::endl
                  << std::endl;

        // set the receive rf gain
        if (vm.count("rx-gain")) {
            std::cout << boost::format("Setting RX Gain: %f dB...") % rx_gain
                      << std::endl;
            // rx_usrp->set_rx_gain(rx_gain, channel);
            rx_usrp->set_normalized_rx_gain(1.0, channel);
            std::cout << boost::format("Actual RX Gain: %f dB...")
                             % rx_usrp->get_rx_gain(channel)
                      << std::endl
                      << std::endl;
        }

        // set the receive analog frontend filter bandwidth
        if (vm.count("rx-bw")) {
            std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (rx_bw / 1e6)
                      << std::endl;
            rx_usrp->set_rx_bandwidth(rx_bw, channel);
            std::cout << boost::format("Actual RX Bandwidth: %f MHz...")
                             % (rx_usrp->get_rx_bandwidth(channel) / 1e6)
                      << std::endl
                      << std::endl;
        }

        // set the receive antenna
        if (vm.count("rx-ant"))
            rx_usrp->set_rx_antenna(rx_ant, channel);
    }

    // Align times in the RX USRP (the TX USRP does not require time-syncing)
    if (rx_usrp->get_num_mboards() > 1) {
        rx_usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
    }

    // for the const wave, set the wave freq for small samples per period
    if (wave_freq == 0 and wave_type == "CONST") {
        wave_freq = tx_usrp->get_tx_rate() / 2;
    }

    // error when the waveform is not possible to generate
    if (std::abs(wave_freq) > tx_usrp->get_tx_rate() / 2) {
        throw std::runtime_error("wave freq out of Nyquist zone");
    }

    // pre-compute the waveform values
    size_t index      = 0;

    // create a transmit streamer
    // linearly map channels (index0 = channel0, index1 = channel1, ...)
    uhd::stream_args_t stream_args("fc32", otw);
    stream_args.channels             = tx_channel_nums;
    uhd::tx_streamer::sptr tx_stream = tx_usrp->get_tx_stream(stream_args);

    // create a receive streamer
    stream_args.channels             = rx_channel_nums;
    uhd::rx_streamer::sptr rx_stream = rx_usrp->get_rx_stream(stream_args);

    // allocate a buffer which we re-use for each channel
    // check samples per buffer, NOTE: default spb is max_spb
    if (spb == 0)
        spb = tx_stream->get_max_num_samps() * 10;
    std::vector<std::complex<float>> buff(spb); 
    int num_channels = tx_channel_nums.size();

    // setup the metadata flags
    uhd::tx_metadata_t tx_md;
    tx_md.start_of_burst = true;
    tx_md.end_of_burst   = false;
    tx_md.has_time_spec  = true;
    tx_md.time_spec = uhd::time_spec_t(0.5); // give us 0.5 seconds to fill the tx buffers

    // Check Ref and LO Lock detect
    std::vector<std::string> tx_sensor_names, rx_sensor_names;
    tx_sensor_names = tx_usrp->get_tx_sensor_names(0);
    if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked")
        != tx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = tx_usrp->get_tx_sensor("lo_locked", 0);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }
    rx_sensor_names = rx_usrp->get_rx_sensor_names(0);
    if (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "lo_locked")
        != rx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = rx_usrp->get_rx_sensor("lo_locked", 0);
        std::cout << boost::format("Checking RX: %s ...") % lo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }

    tx_sensor_names = tx_usrp->get_mboard_sensor_names(0);
    if ((ref == "mimo")
        and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "mimo_locked")
             != tx_sensor_names.end())) {
        uhd::sensor_value_t mimo_locked = tx_usrp->get_mboard_sensor("mimo_locked", 0);
        std::cout << boost::format("Checking TX: %s ...") % mimo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external")
        and (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "ref_locked")
             != tx_sensor_names.end())) {
        uhd::sensor_value_t ref_locked = tx_usrp->get_mboard_sensor("ref_locked", 0);
        std::cout << boost::format("Checking TX: %s ...") % ref_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    rx_sensor_names = rx_usrp->get_mboard_sensor_names(0);
    if ((ref == "mimo")
        and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "mimo_locked")
             != rx_sensor_names.end())) {
        uhd::sensor_value_t mimo_locked = rx_usrp->get_mboard_sensor("mimo_locked", 0);
        std::cout << boost::format("Checking RX: %s ...") % mimo_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(mimo_locked.to_bool());
    }
    if ((ref == "external")
        and (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "ref_locked")
             != rx_sensor_names.end())) {
        uhd::sensor_value_t ref_locked = rx_usrp->get_mboard_sensor("ref_locked", 0);
        std::cout << boost::format("Checking RX: %s ...") % ref_locked.to_pp_string()
                  << std::endl;
        UHD_ASSERT_THROW(ref_locked.to_bool());
    }

    if (total_num_samps == 0) {
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

    // reset usrp time to prepare for transmit/receive
    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    tx_usrp->set_time_now(uhd::time_spec_t(0.0));

    // init rx_buffs
    std::vector<std::vector<std::complex<float>>> rx_buffs(
        rx_channel_nums.size(), std::vector<std::complex<float>>(spb));
    
    // create a vector of pointers to point to each of the channel buffers
    std::vector<std::complex<float>*> buff_ptrs;
    for (size_t i = 0; i < rx_buffs.size(); i++) {
        buff_ptrs.push_back(&rx_buffs[i].front());
    }

    // init tx_buffs
    int seqlen = 137;
    std::vector<std::complex<float>> tx_buff(seqlen*100);
    for (size_t i = 0; i < tx_buff.size(); i++) {
        // compute ZC sequence
        tx_buff[i].real(zc_real[i%seqlen]);
        tx_buff[i].imag(zc_imag[i%seqlen]);
        std::cout << tx_buff[i] << "\n";
    }
    std::vector<std::complex<float>*> tx_buffs(num_channels, &tx_buff.front());
    
    // rx_metadate
    uhd::rx_metadata_t rx_md;

    // init timeout
    double timeout = settling + 0.5f;

    bool overflow_message = true;

    // setup rx_streaming
    uhd::stream_cmd_t stream_cmd((total_num_samps == 0)
                                     ? uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS
                                     : uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    stream_cmd.num_samps  = total_num_samps;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec  = rx_usrp->get_time_now() + uhd::time_spec_t(settling);
    rx_stream->issue_stream_cmd(stream_cmd);

    // on -> true when receive beacon
    timestamp_t t1, t2, t3;

    std::cout << "Samples per buff: " << spb << std::endl;

    tx_stream->send(tx_buffs, seqlen*100, tx_md);
    tx_md.start_of_burst = false;
    tx_md.has_time_spec  = false;

    while (not stop_signal_called)
    {
        /* rx signal */
        size_t num_rx_samps = rx_stream->recv(buff_ptrs, spb, rx_md, timeout);
        t1 = NOW();

        // std::cout << "Num of rcv samples: " << num_rx_samps << std::endl;
        if (std::abs(rx_buffs[0][num_rx_samps-1]) > 0.8) {
            break;
        }
        timeout             = 0.1f; // small timeout for subsequent recv

    }

    t2 = NOW();
    tx_stream->send(tx_buffs, seqlen*100, tx_md);
    t3 = NOW();

    while (not stop_signal_called) {
        /* tx signal */
        // t3 = NOW();
        tx_stream->send(tx_buffs, seqlen*100, tx_md);
    }
    std::cout << "Reaction interval: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() << "ns\n"
            << std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count() << "ns\n";
    // Shut down receiver
    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);

    // transimitter send a mini EOB packet
    tx_md.end_of_burst = true;
    tx_stream->send("", 0, tx_md);

    return EXIT_SUCCESS;
}
