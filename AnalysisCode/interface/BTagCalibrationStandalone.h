#ifndef BTagEntry_H
#define BTagEntry_H

/**
 *
 * BTagEntry
 *
 * Represents one pt- or discriminator-dependent calibration function.
 *
 * measurement_type:    e.g. comb, ttbar, di-mu, boosted, ...
 * sys_type:            e.g. central, plus, minus, plus_JEC, plus_JER, ...
 *
 * Everything is converted into a function, as it is easiest to store it in a
 * txt or json file.
 *
 ************************************************************/

#include <string>
#include <iostream>
#include <exception>
#include <algorithm>
#include <sstream>
#include "TF1.h"
#include "TH1.h"
#include "TSystem.h"


class BTagEntry
{
public:
  enum OperatingPoint {
    OP_LOOSE=0,
    OP_MEDIUM=1,
    OP_TIGHT=2,
    OP_RESHAPING=3,
  };
  enum JetFlavor {
    FLAV_B=0,
    FLAV_C=1,
    FLAV_UDSG=2,
  };
  struct Parameters {
    OperatingPoint operatingPoint;
    std::string measurementType;
    std::string sysType;
    JetFlavor jetFlavor;
    float etaMin;
    float etaMax;
    float ptMin;
    float ptMax;
    float discrMin;
    float discrMax;

    // default constructor
    Parameters(
      OperatingPoint op=OP_TIGHT,
      std::string measurement_type="comb",
      std::string sys_type="central",
      JetFlavor jf=FLAV_B,
      float eta_min=-99999.,
      float eta_max=99999.,
      float pt_min=0.,
      float pt_max=99999.,
      float discr_min=0.,
      float discr_max=99999.
    );

  };

  BTagEntry() {}
  BTagEntry(const std::string &csvLine);
  BTagEntry(const std::string &func, Parameters p);
  BTagEntry(const TF1* func, Parameters p);
  BTagEntry(const TH1* histo, Parameters p);
  ~BTagEntry() {}
  static std::string makeCSVHeader();
  std::string makeCSVLine() const;
  static std::string trimStr(std::string str);

  // public, no getters needed
  std::string formula;
  Parameters params;

};

#endif  // BTagEntry_H


#ifndef BTagCalibration_H
#define BTagCalibration_H

/**
 * BTagCalibration
 *
 * The 'hierarchy' of stored information is this:
 * - by tagger (BTagCalibration)
 *   - by operating point or reshape bin
 *     - by jet parton flavor
 *       - by type of measurement
 *         - by systematic
 *           - by eta bin
 *             - as 1D-function dependent of pt or discriminant
 *
 ************************************************************/

#include <map>
#include <vector>
#include <string>
#include <istream>
#include <ostream>


class BTagCalibration
{
public:
  BTagCalibration() {}
  BTagCalibration(const std::string &tagger);
  BTagCalibration(const std::string &tagger, const std::string &filename);
  ~BTagCalibration() {}

  std::string tagger() const {return tagger_;}

  void addEntry(const BTagEntry &entry);
  const std::vector<BTagEntry>& getEntries(const BTagEntry::Parameters &par) const;

  void readCSV(std::istream &s);
  void readCSV(const std::string &s);
  void makeCSV(std::ostream &s) const;
  std::string makeCSV() const;

protected:
  static std::string token(const BTagEntry::Parameters &par);

  std::string tagger_;
  std::map<std::string, std::vector<BTagEntry> > data_;

};

#endif  // BTagCalibration_H

#ifndef BTagCalibrationReader_H
#define BTagCalibrationReader_H

/**
 * BTagCalibrationReader
 *
 * Helper class to pull out a specific set of BTagEntry's out of a
 * BTagCalibration. TF1 functions are set up at initialization time.
 *
 ************************************************************/

#include <map>
#include <string>
#include <vector>
#include <TF1.h>


class BTagCalibrationReader
{
public:
  BTagCalibrationReader() {}
  BTagCalibrationReader(const BTagCalibration* c,
                        BTagEntry::OperatingPoint op,
                        std::string measurementType="comb",
                        std::string sysType="central");
  ~BTagCalibrationReader() {}

  double eval(BTagEntry::JetFlavor jf,
              float eta,
              float pt,
              float discr=0.) const;

protected:
  struct TmpEntry {
    float etaMin;
    float etaMax;
    float ptMin;
    float ptMax;
    float discrMin;
    float discrMax;
    TF1 func;
  };
  void setupTmpData(const BTagCalibration* c);

  BTagEntry::Parameters params;
  std::map<BTagEntry::JetFlavor, std::vector<TmpEntry> > tmpData_;
  std::vector<bool> useAbsEta;
};

#endif  // BTagCalibrationReader_H
