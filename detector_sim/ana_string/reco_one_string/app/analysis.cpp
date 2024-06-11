#include <fstream>
#include <iostream>
#include <ctime>
#include <sstream>

#include "TTree.h"
#include "TFile.h"

#include "trident/fit.hpp"
#include "trident/filter.hpp"
#include "trident/process_g4_npe.hpp"
#include "trident/root_editor.hpp"
#include "trident/TridentHits.hpp"
#include "Logging.hpp"

// Logging macro to include file, line number, and function name
#define LOG_INFO(message, ...) Logger::getInstance().info("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)
#define LOG_WARN(message, ...) Logger::getInstance().warn("[{}:{}] " message, Logger::getRelativePath(__FILE__), __LINE__, ##__VA_ARGS__)

std::string get_date()
{
    // Get the current time
    std::time_t now = std::time(nullptr);
    // Convert to local time format
    std::tm* local_time = std::localtime(&now);
    // Create a string stream to format the date
    std::ostringstream date_stream;
    date_stream << std::put_time(local_time, "%Y-%m-%d-");
    // Convert to string
    std::string date_str = date_stream.str();

    return date_str;

}
int main(int argc, char** argv){

    // Set the root path to the current directory
    Logger::setRootPath("/Users/meihualin/Projects/trident/siMu_atm/detector_sim/ana_string/reco_one_string"); 

    // Set the log directory before using the logger
    Logger::setLogDirectory("/Users/meihualin/Projects/trident/siMu_atm/detector_sim/ana_string/reco_one_string/logs/");

    int nArgs = 2;
    std::string treename = "Hits";
    std::string tag = "test_looper";
    std::string outfilename = get_date() + tag ;

    // Set a custom logging pattern
    Logger::getInstance().setPattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");

    // Change log levels dynamically
    Logger::getInstance().setConsoleLogLevel(spdlog::level::debug);
    Logger::getInstance().setFileLogLevel(spdlog::level::debug);
    Logger::getInstance().setLoggerLogLevel(spdlog::level::trace);

    if (argc != nArgs) {
        LOG_INFO("Number of argument wrong, expect {}, get {}!", nArgs, argc-1);
        return -1;
    }

    std::string filename = argv[1];

    std::shared_ptr<TFile> root_file = std::make_shared<TFile>(filename.c_str());//, "UPDATE");
    if(root_file->IsZombie()){
        std::cerr << "Error: cannot open file " << "copy.root" << std::endl;
        return -1;
    }

    LOG_INFO("File {} opened.", filename);

    TTree* tree = root_editor::load_tree(root_file.get(), treename);
    LOG_INFO("Tree \"{}\" opened, it contains {} entries.", treename, tree->GetEntries());

    // do something
    LOG_INFO("Outfilename is {}", outfilename);
    TridentHits tridentHits = TridentHits(tree);
    tridentHits.SetOutputname(outfilename);
    tridentHits.Loop(); 

    root_file->Close();

    LOG_INFO("File {} closed.", filename);
    /*
 
    // process_g4_npe npe_process(tree);
    // npe_process.set_qe_file("../qe.csv");
    // npe_process.set_tts(1.3);
    // npe_process.process();
    
    filter filter_obj(tree);
    filter_obj.apply_L2(5);
    TTree* filtered_tree = filter_obj.get_filtered_tree();

    fit::graph_fit(filtered_tree);
   
    root_file->Close();
    */
    return 0;
}                                                                           
