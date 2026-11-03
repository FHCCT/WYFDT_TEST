#pragma once
#include "dt.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

namespace dt {
// CLI replay data. The library's adaptation API and its defaults are unchanged.
inline void readSizeControlProfile(const std::string& path, Args& args, double& globalSize,
    std::unordered_map<int,double>& faceSize, std::unordered_map<int,int>& faceCount,
    std::unordered_map<int,double>& volumeSize, std::unordered_map<int,int>& volumeCount) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot read size-control profile: "+path);
    faceSize.clear(); faceCount.clear(); volumeSize.clear(); volumeCount.clear();
    std::map<std::string,double*> real={{"global_size",&globalSize},{"optTh",&args.optTh},
        {"optanglestrict",&args.optanglestrict},{"adpangle",&args.adpangle}};
    std::map<std::string,int*> integer={{"optloop",&args.optloop},{"optlevel",&args.optlevel}};
    std::set<std::string> seen;
    std::string text; int line=0;
    auto error=[&]() { return std::runtime_error("Invalid size-control profile at "+path+":"+std::to_string(line)); };
    while (std::getline(file,text)) {
        ++line; text=text.substr(0,text.find('#'));
        std::istringstream input(text); std::string key,extra;
        if (!(input>>key)) continue;
        if (real.count(key)) {
            double value;
            if (!(input>>value)||!std::isfinite(value)||!seen.insert(key).second) throw error();
            *real[key]=value;
        } else if (integer.count(key)) {
            int value;
            if (!(input>>value)||!seen.insert(key).second) throw error();
            *integer[key]=value;
        } else if (key=="face_size" || key=="volume_size") {
            int id; double size;
            if (!(input>>id>>size)||id<0||!(size>0)||!std::isfinite(size)) throw error();
            auto& sizes=key=="face_size"?faceSize:volumeSize;
            if (!sizes.emplace(id,size).second) throw error();
        } else if (key=="face_count" || key=="volume_count") {
            int id,count;
            if (!(input>>id>>count)||id<0||count<=0) throw error();
            auto& counts=key=="face_count"?faceCount:volumeCount;
            if (!counts.emplace(id,count).second) throw error();
        } else throw error();
        if (input>>extra) throw error();
    }
    if (!seen.count("global_size") || !(globalSize>0) || args.optloop<1 || args.optlevel<0 ||
        args.optTh<0 || args.optTh>1 || args.optanglestrict<0 || args.optanglestrict>180 ||
        args.adpangle<0 || args.adpangle>180) throw error();
}
}
