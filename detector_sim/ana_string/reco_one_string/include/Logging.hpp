#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <memory>
#include <string>
#include <filesystem>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

class Logger {
public:
    // Get the singleton instance of the Logger class
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }
    
    // Method to log debugings with formatting
    template<typename... Args>
    void debug(const std::string& fmt, Args&&... args) {
        logger->debug(fmt, std::forward<Args>(args)...);
    }

    // Method to log warnings with formatting
    template<typename... Args>
    void warn(const std::string& fmt, Args&&... args) {
        logger->warn(fmt, std::forward<Args>(args)...);
    }

    // Method to log info with formatting
    template<typename... Args>
    void info(const std::string& fmt, Args&&... args) {
        logger->info(fmt, std::forward<Args>(args)...);
    }

    // Method to set the log directory
    static void setLogDirectory(const std::string& directory) {
        logDirectory = directory;
        getInstance().initializeLogger();
    }

    // Method to set the root path for computing relative file paths
    static void setRootPath(const std::string& path) {
        rootPath = std::filesystem::absolute(path);
    }

    // Methods to set log levels dynamically
    void setConsoleLogLevel(spdlog::level::level_enum level) {
        console_sink->set_level(level);
    }

    void setFileLogLevel(spdlog::level::level_enum level) {
        file_sink->set_level(level);
    }

    void setLoggerLogLevel(spdlog::level::level_enum level) {
        logger->set_level(level);
    }

    // Method to set logging pattern
    void setPattern(const std::string& pattern) {
        logger->set_pattern(pattern);
    }

    // Method to get the relative file path
    static std::string getRelativePath(const std::string& filepath) {
        std::filesystem::path absPath = std::filesystem::absolute(filepath);
        return std::filesystem::relative(absPath, rootPath).string();
    }

private:
    std::shared_ptr<spdlog::logger> logger;
    static std::string logDirectory;
    static std::filesystem::path rootPath;
    std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> console_sink;
    std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;

    // Private constructor to prevent instantiation
    Logger() {
        initializeLogger();
    }

    void initializeLogger() {
        console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::warn);
        console_sink->set_pattern("[%^%l%$] %v");

        std::string logFilePath = logDirectory + "/multisink.txt";
        file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFilePath, true);
        file_sink->set_level(spdlog::level::trace);

        logger = std::make_shared<spdlog::logger>("logger", spdlog::sinks_init_list{console_sink, file_sink});
        logger->set_level(spdlog::level::debug);

        // Default pattern, excluding the logger name and using the relative filename
        logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%s:%#] %v");
    }

    // Delete copy constructor and assignment operator
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
};

#endif // LOGGING_HPP

