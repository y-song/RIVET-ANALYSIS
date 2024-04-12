// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"

#include "cmath"
#include "iostream"
#include "fstream"
#include "random"
#include "vector"
#include "stdio.h"
#include "stdlib.h"

using namespace std;

namespace fastjet
{
  class MyUserInfo : public PseudoJet::UserInfoBase
  {
  public:
    MyUserInfo(vector<double> q, const int &matching, const int &parton) : _q(q), _parton(parton) {}
    vector<double> q() const { return _q; }
    int parton() const { return _parton; }

  protected:
    vector<double> _q;
    int _parton;
  };

  class DihadronInfo : public PseudoJet::UserInfoBase
  {
  public:
    DihadronInfo(const double &b_ch, const double &b, const double &y1_ch, const double &y2_ch, const double &phi1_ch, const double &phi2_ch, const double &y1, const double &y2, const double &phi1, const double &phi2, const double &q, const int &parton) : _b_ch(b_ch), _b(b), _y1_ch(y1_ch), _y2_ch(y2_ch), _phi1_ch(phi1_ch), _phi2_ch(phi2_ch), _y1(y1), _y2(y2), _phi1(phi1), _phi2(phi2), _q(q), _parton(parton) {}
    double b_ch() const { return _b_ch; }
    double b() const { return _b; }
    double y1_ch() const { return _y1_ch; }
    double y2_ch() const { return _y2_ch; }
    double phi1_ch() const { return _phi1_ch; }
    double phi2_ch() const { return _phi2_ch; }
    double y1() const { return _y1; }
    double y2() const { return _y2; }
    double phi1() const { return _phi1; }
    double phi2() const { return _phi2; }
    double q() const { return _q; }
    int parton() const { return _parton; }

  protected:
    double _b_ch;
    double _b;
    double _y1_ch;
    double _y2_ch;
    double _phi1_ch;
    double _phi2_ch;
    double _y1;
    double _y2;
    double _phi1;
    double _phi2;
    double _q;
    int _parton;
  };

  class ParticleInfo : public PseudoJet::UserInfoBase
  {
  public:
    ParticleInfo(const double charge, const int &pid) : _charge(charge), _pid(pid) {}
    double charge() const { return _charge; }
    int pid() const { return _pid; }

  protected:
    double _charge;
    int _pid;
  };
}

namespace Rivet
{
  class HERWIG_JETS_rc : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERWIG_JETS_rc);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // load tracking efficiency file and save values to an array
      efftxtfile.open("eff.txt");
      if (!efftxtfile.good())
      {
        fputs("efficiency file not found\n", stderr);
        abort();
      }
      string line;
      for (int i = 0; getline(efftxtfile, line); ++i)
      {
        istringstream iss(line);
        string delimiter = ", ";
        string token;
        size_t pos = 0;
        int j = 0;
        while ((pos = line.find(delimiter)) != string::npos)
        {
          token = line.substr(0, pos);
          eff_array[i][j] = stod(token);
          line.erase(0, pos + delimiter.length());
          j += 1;
        }
      }

      // open output file
      mytxtfile.open("herwig_50pthat_rc_r06.txt");

      // make JP grids
      for (int i = 0; i <= Nbounds_JP_eta; ++i)
      {
        double etax = -1 * etaMax + i * (double)2 * etaMax / Nbounds_JP_eta;
        etabins_JP.push_back(etax);
      }
      for (int i = 0; i <= Nbounds_JP_phi; ++i)
      {
        double phix = i * (double)phiMax / Nbounds_JP_phi;
        phibins_JP.push_back(phix);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {
      int unmatched = -1;
      int matched = 1;
      int missing = 2;
      int fake = 3;
      int isJP2 = 0;
      double JP2_et = -1;
      int JP2_eta = -1;
      int JP2_phi = -1;

      e_JP.clear();
      for (int i = 0; i <= Nbounds_JP_eta; ++i)
      {
        vector<double> xgrid;
        for (int j = 0; j <= Nbounds_JP_phi; ++j)
        {
          xgrid.push_back(0.0);
        }
        e_JP.push_back(xgrid);
        xgrid.clear();
      }

      vector<double> k{2};
      double xsecweight = handler().nominalCrossSection();

      //===========================================
      // process particles
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();
      PseudoJets parts, parts_s;

      for (const Particle &p : fsParticles)
      {
        if (p.pid() == 12 || p.pid() == 14 || p.pid() == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        int pid = p.pid();
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid));
        parts.push_back(pseudojet);

        random_device generator;
        bool neutral = !p.isCharged();
        if (neutral)
        {
          double es, ratio;
          normal_distribution<double> smearNeutral(p.E(), 0.14 * p.E());
          es = smearNeutral(generator);
          ratio = es / p.E();
          PseudoJet pj_smeared(p.px() * ratio, p.py() * ratio, p.pz() * ratio, es);
          pj_smeared.reset_momentum(pj_smeared.px(), pj_smeared.py(), pj_smeared.pz(), sqrt(pow(pj_smeared.px(), 2) + pow(pj_smeared.py(), 2) + pow(pj_smeared.pz(), 2)));
          if (pj_smeared.Et() < 0.2 || pj_smeared.Et() > 30. || abs(pj_smeared.eta()) > 1.0)
            continue;
          pj_smeared.set_user_info(new fastjet::ParticleInfo(charge, pid));
          parts_s.push_back(pj_smeared);
        }
        else
        {
          uniform_real_distribution<double> dropCharge(0.0, 1.0);
          double keep = dropCharge(generator);
          double eff;
          if (p.pt() < 2)
          {
            eff = eff_array[(int)floor((p.pt() - 0.2) / 0.04)][(int)floor((p.eta() + 1.975) / 0.049375)];
          }
          else
          {
            eff = eff_array[44][(int)floor((p.eta() + 1.975) / 0.049375)];
          }
          if (keep > eff)
            continue;
          double pts, ratio;
          normal_distribution<double> smearCharge(p.pt(), -0.026 + 0.02 * p.pt() + 0.003 * pow(p.pt(), 2));
          pts = smearCharge(generator);
          ratio = pts / p.pt();
          PseudoJet pj_smeared(p.px() * ratio, p.py() * ratio, p.pz() * ratio, p.E() * ratio);
          pj_smeared.reset_momentum(pj_smeared.px(), pj_smeared.py(), pj_smeared.pz(), sqrt(pow(pj_smeared.px(), 2) + pow(pj_smeared.py(), 2) + pow(pj_smeared.pz(), 2) + pow(0.1395704, 2))); // charged pion mass
          if (pj_smeared.perp() < 0.2 || pj_smeared.perp() > 30. || abs(pj_smeared.eta()) > 1.0)
            continue;
          pj_smeared.set_user_info(new fastjet::ParticleInfo(charge, pid));
          parts_s.push_back(pj_smeared);
        }
      }

      //===========================================
      // check if event passes JP2 trigger
      for (PseudoJet &pj_smeared : parts_s)
      {
        for (int x = 0; x < Nbounds_JP_eta; ++x)
        {
          for (int y = 0; y < Nbounds_JP_phi; ++y)
          {
            // neutrals or electrons
            if (pj_smeared.user_info<fastjet::ParticleInfo>().charge() == 0 || abs(pj_smeared.user_info<fastjet::ParticleInfo>().pid() == 11))
            {
              if (pj_smeared.eta() > etabins_JP[x] && pj_smeared.eta() < etabins_JP[x + 1] && pj_smeared.phi() > phibins_JP[y] && pj_smeared.phi() < phibins_JP[y + 1])
              {
                e_JP[x][y] += pj_smeared.Et();
                break;
              }
            }
            // charged particles with pT > 0.35 GeV
            else if (pj_smeared.perp() > 0.35)
            {
              if (pj_smeared.eta() > etabins_JP[x] && pj_smeared.eta() < etabins_JP[x + 1] && pj_smeared.phi() > phibins_JP[y] && pj_smeared.phi() < phibins_JP[y + 1])
              {
                e_JP[x][y] += 0.35;
                break;
              }
            }
          }
        }
      }

      for (int x = 0; x < Nbounds_JP_eta; ++x)
      {
        if (isJP2 == 1)
          break;
        for (int y = 0; y < Nbounds_JP_phi; ++y)
        {
          if (e_JP[x][y] > 7.3)
          {
            isJP2 = 1;
            JP2_et = e_JP[x][y];
            JP2_eta = x;
            JP2_phi = y;
            break;
          }
        }
      }

      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      // int evid = theEvent->event_number();

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(5.0) * fastjet::SelectorEtaRange(-0.4, 0.4);
      fastjet::Selector selector_s = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.4, 0.4); // * fastjet::SelectorMassMin(1.0);

      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.6);

      fastjet::contrib::SoftDrop sd(0, 0.1);

      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
      vector<fastjet::PseudoJet> sdjets = sd(jets);
      // vector<fastjet::PseudoJet> partons = {};

      fastjet::ClusterSequence cs_s(parts_s, jet_def);
      vector<fastjet::PseudoJet> jets_s = sorted_by_pt(selector_s(cs_s.inclusive_jets()));
      vector<fastjet::PseudoJet> sdjets_s = sd(jets_s);

      //-------------------------------------------
      // calculate pythia-level jet charge and leading dihadrons
      double kappa = 0.5;
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        double b_ch = 0;
        double b = 0;
        double y1_ch = 0;
        double y2_ch = 0;
        double phi1_ch = 0;
        double phi2_ch = 0;
        double y1 = 0;
        double y2 = 0;
        double phi1 = 0;
        double phi2 = 0;
        double q = 0;
        vector<fastjet::PseudoJet> IncPart = sorted_by_pt(jets[i].constituents());
        vector<fastjet::PseudoJet> ChargedPart;
        if (IncPart.size() < 2)
        {
          b = -9;
          y1 = -9;
          y2 = -9;
          phi1 = -9;
          phi2 = -9;
        }
        else
        {
          b = (IncPart.at(0)).user_info<fastjet::ParticleInfo>().charge() * (IncPart.at(1)).user_info<fastjet::ParticleInfo>().charge();
          y1 = (IncPart.at(0)).rap();
          y2 = (IncPart.at(1)).rap();
          phi1 = (IncPart.at(0)).phi();
          phi2 = (IncPart.at(1)).phi();
        }
       double numerator = 0;
        for (unsigned int j = 0; j < IncPart.size(); j++)
        {
          PseudoJet part = IncPart.at(j);
          double charge = part.user_info<fastjet::ParticleInfo>().charge();
          if (charge != 0)
          {
            part.set_user_info(new fastjet::ParticleInfo(charge, -999));
            ChargedPart.push_back(part);
          }
          numerator += pow(part.perp(), kappa) * charge;
        }
        q = numerator / pow(jets[i].perp(), kappa);

        if (ChargedPart.size() < 2)
        {
          b_ch = -9;
          y1_ch = -9;
          y2_ch = -9;
          phi1_ch = -9;
          phi2_ch = -9;
        }
        else
        {
          b_ch = (ChargedPart.at(0)).user_info<fastjet::ParticleInfo>().charge() * (ChargedPart.at(1)).user_info<fastjet::ParticleInfo>().charge();
          y1_ch = (ChargedPart.at(0)).rap();
          y2_ch = (ChargedPart.at(1)).rap();
          phi1_ch = (ChargedPart.at(0)).phi();
          phi2_ch = (ChargedPart.at(1)).phi();
        }
        jets[i].set_user_index(-1);
        jets[i].set_user_info(new fastjet::DihadronInfo(b_ch, b, y1_ch, y2_ch, phi1_ch, phi2_ch, y1, y2, phi1, phi2, q, -999));
        IncPart.clear();
        ChargedPart.clear();
      }
      
      //-------------------------------------------
      // calculate smeared jet charge and dihadrons
      for (unsigned int is = 0; is < jets_s.size(); is++)
      {
        double bs_ch = 0;
        double bs = 0;
        double y1s_ch = 0;
        double y2s_ch = 0;
        double phi1s_ch = 0;
        double phi2s_ch = 0;
        double y1s = 0;
        double y2s = 0;
        double phi1s = 0;
        double phi2s = 0;
        double qs = 0;
 	vector<fastjet::PseudoJet> IncPart = sorted_by_pt(jets_s[is].constituents());
        vector<fastjet::PseudoJet> ChargedPart;
        if (IncPart.size() < 2)
        {
          bs = -9;
          y1s = -9;
          y2s = -9;
          phi1s = -9;
          phi2s = -9;
        }
        else
        {
          bs = (IncPart.at(0)).user_info<fastjet::ParticleInfo>().charge() * (IncPart.at(1)).user_info<fastjet::ParticleInfo>().charge();
          y1s = (IncPart.at(0)).rap();
          y2s = (IncPart.at(1)).rap();
          phi1s = (IncPart.at(0)).phi();
          phi2s = (IncPart.at(1)).phi();
        }
       double numerator = 0;
        for (unsigned int j = 0; j < IncPart.size(); j++)
        {
          PseudoJet part = IncPart.at(j);
          double charge = part.user_info<fastjet::ParticleInfo>().charge();
          if (charge != 0)
          {
            part.set_user_info(new fastjet::ParticleInfo(charge, -999));
            ChargedPart.push_back(part);
          }
          numerator += pow(part.perp(), kappa) * charge;
        }
        qs = numerator / pow(jets_s[is].perp(), kappa);
      
      if (ChargedPart.size() < 2)
      {
        bs_ch = -9;
        y1s_ch = -9;
        y2s_ch = -9;
        phi1s_ch = -9;
        phi2s_ch = -9;
      }
      else
      {
        bs_ch = (ChargedPart.at(0)).user_info<fastjet::ParticleInfo>().charge() * (ChargedPart.at(1)).user_info<fastjet::ParticleInfo>().charge();
        y1s_ch = (ChargedPart.at(0)).rap();
        y2s_ch = (ChargedPart.at(1)).rap();
        phi1s_ch = (ChargedPart.at(0)).phi();
        phi2s_ch = (ChargedPart.at(1)).phi();
      }
      jets_s[is].set_user_index(-1);
      jets_s[is].set_user_info(new fastjet::DihadronInfo(bs_ch, bs, y1s_ch, y2s_ch, phi1s_ch, phi2s_ch, y1s, y2s, phi1s, phi2s, qs, -999));
      IncPart.clear();
      ChargedPart.clear();
    }

    //-------------------------------------------
    // Match smeared jets with unsmeared jets
    for (unsigned int i = 0; i < jets.size(); i++)
    {
      for (unsigned int is = 0; is < jets_s.size(); is++)
      {
        if (jets_s[is].user_index() == matched)
          continue;
        if (jets[i].delta_R(jets_s[is]) < 0.4)
        {

          jets[i].set_user_index(matched);
          jets_s[is].set_user_index(matched);
          mytxtfile << jets[i].user_info<fastjet::DihadronInfo>().b_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().b() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y2() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi2() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().b_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().b() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y1_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y2_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi1_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi2_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y1() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y2() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi1() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi2() << ", " << jets[i].user_info<fastjet::DihadronInfo>().q() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().q() << ", ";

          /*for (vector<double>::iterator it = jet_q.begin(); it != jet_q.end(); ++it)
          {
            mytxtfile << *it << ", ";
          }
          for (vector<double>::iterator it = jet_q_s.begin(); it != jet_q_s.end(); ++it)
          {
            mytxtfile << *it << ", ";
          }*/

          mytxtfile << xsecweight << ", " << matched << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << -999 << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << jets_s[is].perp() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << -999 << ", " << sdjets_s[is].m() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << i << ", " << theEvent->event_scale() << "\n";
          break;
        }
      }
    }

    for (unsigned int i = 0; i < jets.size(); i++)
    {
      if (jets[i].user_index() == matched)
      {
        continue;
      }
      jets[i].set_user_index(missing);

      mytxtfile << jets[i].user_info<fastjet::DihadronInfo>().b_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().b() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y2() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi2() << ", -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, " << jets[i].user_info<fastjet::DihadronInfo>().q() << ", -9, ";

      /*for (vector<double>::iterator it = jet_q.begin(); it != jet_q.end(); ++it)
      {
        mytxtfile << *it << ", -999, ";
      }*/

      mytxtfile << xsecweight << ", " << missing << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << -999 << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << i << ", " << theEvent->event_scale() << "\n";
    }
    for (unsigned int is = 0; is < jets_s.size(); is++)
    {
      if (jets_s[is].user_index() == matched)
        continue;
      jets_s[is].set_user_index(fake);

      mytxtfile << "-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, " << jets_s[is].user_info<fastjet::DihadronInfo>().b_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().b() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y1_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y2_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi1_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi2_ch() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y1() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().y2() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi1() << ", " << jets_s[is].user_info<fastjet::DihadronInfo>().phi2() << ", -9, " << jets_s[is].user_info<fastjet::DihadronInfo>().q() << ", ";

      /*mytxtfile << "-999, ";
      for (vector<double>::iterator it = jet_q_s.begin(); it != jet_q_s.end(); ++it)
      {
        mytxtfile << *it << ", ";
      }*/

      mytxtfile << xsecweight << ", " << fake << ", -999, -999, -999, -999, -999, -999, -999, " << jets_s[is].perp() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", -999, " << sdjets_s[is].m() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", -999, -999\n";
    }
  }

  /// Normalise histograms etc., after the run
  void finalize()
  {
    double norm = crossSection() / sumOfWeights();
    std::cout << norm << std::endl;
    mytxtfile.close();
  }

private:
  double eff_array[45][80];
  std::ofstream mytxtfile;
  ifstream efftxtfile;
  double etaMax = 1.0;
  double pi = 3.14159265359;
  double phiMax = 2 * 3.14159265359;
  int Nbounds_JP_eta = (int)2 * etaMax;
  int Nbounds_JP_phi = (int)phiMax;
  vector<double> etabins_JP;
  vector<double> phibins_JP;
  vector<vector<double>> e_JP;
};

RIVET_DECLARE_PLUGIN(HERWIG_JETS_rc);
}
