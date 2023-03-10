/*
 * ae_logger.h
 *
 *  Created on: Jan 26, 2015
 *      Author: arrouan
 */

#ifndef AE_LOGGER_H_
#define AE_LOGGER_H_

#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <vector>
#include <sstream>
#include <mutex>

using namespace std;

enum logger_category {
	SELECTION = 0,
	TOTAL = 1,
  SEARCH_START_RNA = 2,
  DNA_TO_RNA = 3,
  RNA_TO_PROTEIN = 4,
  COMPUTE_PROTEIN = 5,
  PROTEIN_TO_PHENOTYPE = 6
};

class ae_logger
{
  public :
	static inline void init(string file);
	static inline void addLog(int category, string msg);
	static inline void addLog(int category, long int msg);
	static inline void flush(int generation);
	static unordered_map<int,unordered_multiset<string>> logMap;
	static string logFile;
	static std::mutex loggerMtx;
};

//int cpt = 0;

void ae_logger::init(string file) {
	ae_logger::logFile = file;
}

void ae_logger::addLog(int category, string msg) {
	loggerMtx.lock();
	ae_logger::logMap[category].insert(msg);
	loggerMtx.unlock();
}

void ae_logger::addLog(int category, long int msg) {
	std::stringstream ss;
	ss << msg;
	string msg_str = ss.str();
	loggerMtx.lock();
	ae_logger::logMap[category].insert(msg_str);
	loggerMtx.unlock();
}

void ae_logger::flush(int generation) {
	ofstream loggerFile;
	loggerFile.open(ae_logger::logFile,ios::out | ios::app );

//	printf("CPT: %d\n",cpt);cpt=0;
	loggerMtx.lock();
  long int avg = 0;
  int cpt = 0;
	for (int i = 0; i <= 28; i++) {
    avg = 0;
    cpt = 0;
    switch (i) {
      case SELECTION:
        loggerFile << "SELECTION," << generation;
        for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
          loggerFile << "," << *it;
        }
        loggerFile << endl;
        logMap[i].clear();
        break;
      case TOTAL:
        if (logMap[i].size() > 0) {
          loggerFile << "TOTAL," << generation;
          for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
            loggerFile << "," << *it;
          }
          loggerFile << endl;
          logMap[i].clear();
        }
        break;
      case SEARCH_START_RNA:
        if (logMap[i].size() > 0) {
          for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
            avg += std::stol(*it);
            cpt++;
          }
          avg/=cpt;
          loggerFile<<"SEARCH_START_RNA,"<<generation<<","<<avg<<endl;
          logMap[i].clear();
        }
        break;
      case DNA_TO_RNA:
        if (logMap[i].size() > 0) {
          for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
            avg += std::stol(*it);
            cpt++;
          }
          avg/=cpt;
          loggerFile<<"DNA_TO_RNA,"<<generation<<","<<avg<<endl;
          logMap[i].clear();
        }
        break;
      case RNA_TO_PROTEIN:
        if (logMap[i].size() > 0) {
          for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
            avg += std::stol(*it);
            cpt++;
          }
          avg/=cpt;
          loggerFile<<"RNA_TO_PROTEIN,"<<generation<<","<<avg<<endl;
          logMap[i].clear();
        }
        break;
      case COMPUTE_PROTEIN:
        if (logMap[i].size() > 0) {
          for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
            avg += std::stol(*it);
            cpt++;
          }
          avg/=cpt;
          loggerFile<<"COMPUTE_PROTEIN,"<<generation<<","<<avg<<endl;
          logMap[i].clear();
        }
        break;
      case PROTEIN_TO_PHENOTYPE:
        if (logMap[i].size() > 0) {
          for (auto it = logMap[i].begin(); it != logMap[i].end(); ++it) {
            avg += std::stol(*it);
            cpt++;
          }
          avg/=cpt;
          loggerFile<<"PROTEIN_TO_PHENOTYPE,"<<generation<<","<<avg<<endl;
          logMap[i].clear();
        }
        break;
    }
  }
  loggerMtx.unlock();
}

#endif /* AE_LOGGER_H_ */
