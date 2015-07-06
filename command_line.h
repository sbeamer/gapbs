// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef COMMAND_LINE_H_
#define COMMAND_LINE_H_

#include <getopt.h>
#include <stdio.h>

#include <algorithm>
#include <cstdlib>
#include <cinttypes>
#include <iostream>
#include <string>
#include <vector>


class CLBase {
 protected:
  int argc_;
  char** argv_;
  std::string name_;
  std::string get_args_ = "f:g:hsu:";
  std::vector<std::string> help_strings_;

  bool ok_to_continue_ = true;
  int scale_ = -1;
  std::string filename_ = "";
  bool symmetrize_ = false;
  bool uniform_ = false;

  void AddHelpLine(char opt, std::string opt_arg, std::string text,
                   std::string def="") {
    const int buf_len = 100;
    char buf[buf_len];
    if (opt_arg != "")
      opt_arg = "<" + opt_arg + ">";
    if (def != "")
      def = "[" + def + "]";
    sprintf(buf, " -%c %-9s: %-57s%7s", opt, opt_arg.c_str(),
            text.c_str(), def.c_str());
    help_strings_.push_back(buf);
  }

 public:
  CLBase(int argc, char** argv, std::string name="") :
         argc_(argc), argv_(argv), name_(name) {
    AddHelpLine('h', "", "print this help message");
    AddHelpLine('f', "file", "load graph from file");
    AddHelpLine('s', "", "symmetrize input edge list", "false");
    AddHelpLine('g', "scale", "generate 2^scale kronecker graph");
    AddHelpLine('u', "scale", "generate 2^scale uniform-random graph");
  }

  bool ParseArgs() {
    signed char c_opt;
    extern char *optarg;          // from and for getopt
    while((c_opt = getopt(argc_, argv_, get_args_.c_str())) != -1) {
      HandleArg(c_opt, optarg);
    }
    if ((filename_ == "") && (scale_ == -1)) {
      std::cout << "No graph input specified. (Use -h for help)" << std::endl;
      ok_to_continue_ = false;
    }
    if (scale_ != -1)
      symmetrize_ = true;
    return ok_to_continue_;
  }

  void virtual HandleArg(signed char opt, char* opt_arg) {
    switch (opt) {
      case 'f': filename_ = std::string(opt_arg);           break;
      case 'g': scale_ = atoi(opt_arg);                     break;
      case 'h': PrintUsage();                               break;
      case 's': symmetrize_ = true;                         break;
      case 'u': uniform_ = true; scale_ = atoi(opt_arg);    break;
    }
  }

  void PrintUsage() {
    std::cout << name_ << std::endl;
    // std::sort(help_strings_.begin(), help_strings_.end());
    for (std::string h : help_strings_)
      std::cout << h << std::endl;
    ok_to_continue_ = false;
  }

  int scale() { return scale_; }
  std::string filename() { return filename_; }
  bool symmetrize() { return symmetrize_; }
  bool uniform() { return uniform_; }
};



class CLApp : public CLBase {
  bool do_analysis_ = false;
  int num_trials_ = 16;
  int64_t start_vertex_ = -1;

 public:
  CLApp(int argc, char** argv, std::string name) : CLBase(argc, argv, name) {
    get_args_ += "an:r:";
    char buf[30];
    sprintf(buf, "%d", num_trials_);
    AddHelpLine('a', "a", "output analysis of last run", "false");
    AddHelpLine('n', "n", "perform n trials", buf);
    AddHelpLine('r', "node", "start from node r", "rand");
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'a': do_analysis_ = true;                    break;
      case 'n': num_trials_ = atoi(opt_arg);            break;
      case 'r': start_vertex_ = atol(opt_arg);          break;
      default: CLBase::HandleArg(opt, opt_arg);
    }
  }

  bool do_analysis() { return do_analysis_; }
  int num_trials() { return num_trials_; }
  int64_t start_vertex() { return start_vertex_; }
};



class CLIterApp : public CLApp {
  int num_iters_;

 public:
  CLIterApp(int argc, char** argv, std::string name, int num_iters) :
    CLApp(argc, argv, name), num_iters_(num_iters) {
    get_args_ += "k:";
    char buf[15];
    sprintf(buf, "%d", num_iters_);
    AddHelpLine('k', "k", "perform k iterations", buf);
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'k': num_iters_ = atoi(opt_arg);            break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  int num_iters() { return num_iters_; }
};



class CLDelta : public CLApp {
  int delta_ = 1;

 public:
  CLDelta(int argc, char** argv, std::string name) : CLApp(argc, argv, name) {
    get_args_ += "d:";
    char buf[15];
    sprintf(buf, "%d", delta_);
    AddHelpLine('d', "d", "delta parameter", buf);
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'd': delta_ = atoi(opt_arg);               break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  int delta() { return delta_; }
};



class CLConvert : public CLBase {
  std::string out_filename_ = "";
  bool out_weighted_ = false;
  bool out_el_ = false;
  bool out_sg_ = false;

 public:
  CLConvert(int argc, char** argv, std::string name)
      : CLBase(argc, argv, name) {
    get_args_ += "e:b:w";
    AddHelpLine('b', "file", "output serialized graph to file");
    AddHelpLine('e', "file", "output edge list to file");
    AddHelpLine('w', "file", "make output weighted");
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'b': out_sg_=true; out_filename_=std::string(opt_arg);   break;
      case 'e': out_el_=true; out_filename_=std::string(opt_arg);   break;
      case 'w': out_weighted_ = true;                               break;
      default: CLBase::HandleArg(opt, opt_arg);
    }
  }

  std::string out_filename() { return out_filename_; }
  bool out_weighted() { return out_weighted_; }
  bool out_el() { return out_el_; }
  bool out_sg() { return out_sg_; }
};


#endif  // COMMAND_LINE_H_
