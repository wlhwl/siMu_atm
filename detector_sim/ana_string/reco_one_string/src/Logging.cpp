#include "Logging.hpp"

// Initialize the static member
std::string Logger::logDirectory = "logs";
std::filesystem::path Logger::rootPath = std::filesystem::current_path();

