#pragma once
#include <codecvt>
#include "3rd/SomeLogger.h"
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#define M_PI 3.14159265358979323846

#define E_LOG(s) SomeLogger::Logger::Instance().printEndl(true).log(SomeLogger::LoggerLevel::ERR, SomeLogger::Color::Red, SomeLogger::Color::Black) << s;
#define I_LOG(s) SomeLogger::Logger::Instance().printEndl(true).log(SomeLogger::LoggerLevel::INFO, SomeLogger::Color::Green, SomeLogger::Color::Black) << s;
#define D_LOG(s) SomeLogger::Logger::Instance().printEndl(true).log(SomeLogger::LoggerLevel::DEBUG, SomeLogger::Color::Debug, SomeLogger::Color::Black) << s;
#define W_LOG(s) SomeLogger::Logger::Instance().printEndl(true).log(SomeLogger::LoggerLevel::WARNING, SomeLogger::Color::Red, SomeLogger::Color::Black) << s;

std::string ws2s(const std::wstring& wstr) {
    using convert_typeX = std::codecvt_utf8<wchar_t>;
    std::wstring_convert<convert_typeX, wchar_t> converterX;

    return converterX.to_bytes(wstr);
}

struct CenterPoint {
    double x = 0.0, y = 0.0;
    CenterPoint(double _x=0.0, double _y=0.0): x(_x), y(_y) {}
    ~CenterPoint() = default;
};

std::vector<std::string> split(const std::string& str, char delim = ' ') {
    std::stringstream ss(str);
    std::string token;
    std::vector<std::string> cont;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
	return cont;
}
