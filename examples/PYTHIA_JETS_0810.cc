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
    MyUserInfo(vector<double> q, const int &matching, const int &parton) : _q(q), _matching(matching), _parton(parton) {}
    vector<double> q() const { return _q; }
    int matching() const { return _matching; }
    int parton() const { return _parton; }

  protected:
    vector<double> _q;
    int _matching;
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
  class PYTHIA_JETS_0810 : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_0810);

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
      mytxtfile.open("pythia1208_50pthat.txt");

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

      vector<double> k{0.5, 2};
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
      int evid = theEvent->event_number();

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(5.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::Selector selector_s = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.6, 0.6) * fastjet::SelectorMassMin(1.0);

      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);

      fastjet::contrib::SoftDrop sd(0, 0.1);

      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
      vector<fastjet::PseudoJet> sdjets = sd(jets);
      vector<fastjet::PseudoJet> partons = {};

      fastjet::ClusterSequence cs_s(parts_s, jet_def);
      vector<fastjet::PseudoJet> jets_s = sorted_by_pt(selector_s(cs_s.inclusive_jets()));
      vector<fastjet::PseudoJet> sdjets_s = sd(jets_s);

      //-------------------------------------------
      // calculate pythia-level jet charge
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        vector<double> q{};
        for (const double &kappa : k)
        {
          double numerator = 0;
          for (unsigned int j = 0; j < jets[i].constituents().size(); j++)
          {
            PseudoJet part = jets[i].constituents().at(j);
            double charge = part.user_info<fastjet::ParticleInfo>().charge();
            numerator += pow(part.perp(), kappa) * charge;
          }
          q.push_back(numerator / pow(jets[i].perp(), kappa));
        }
        jets[i].set_user_info(new fastjet::MyUserInfo(q, -1, -999));
      }

      //-------------------------------------------
      // calculate smeared jet charge
      for (unsigned int is = 0; is < jets_s.size(); is++)
      {
        vector<double> q{};
        for (const double &kappa : k)
        {
          double numerator = 0;
          for (unsigned int j = 0; j < jets_s[is].constituents().size(); j++)
          {
            PseudoJet part = jets_s[is].constituents().at(j);
            double charge = part.user_info<fastjet::ParticleInfo>().charge();
            numerator += pow(part.perp(), kappa) * charge;
          }
          q.push_back(numerator / pow(jets_s[is].perp(), kappa));
        }
        jets_s[is].set_user_info(new fastjet::MyUserInfo(q, -1, -999));
      }

     //-------------------------------------------
     // Match jets to partons
      for (HepMC::GenEvent::particle_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
      {

        if ((*p)->status() != 23)
          continue;
        PseudoJet parton_pseudojet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        partons.push_back(parton_pseudojet);

        for (unsigned int i = 0; i < jets.size(); i++)
        {
          if (jets[i].delta_R(parton_pseudojet) > 0.4 || jets[i].user_info<fastjet::MyUserInfo>().parton() != -999)
            continue;
          vector<double> jet_q = jets[i].user_info<fastjet::MyUserInfo>().q();
          int matching = jets[i].user_info<fastjet::MyUserInfo>().matching();
          jets[i].set_user_info(new fastjet::MyUserInfo(jet_q, matching, (*p)->pdg_id()));
          break;
        }

        for (unsigned int is = 0; is < jets_s.size(); is++)
        {
          if (jets_s[is].delta_R(parton_pseudojet) > 0.4 || jets_s[is].user_info<fastjet::MyUserInfo>().parton() != -999)
            continue;
          vector<double> jet_q_s = jets_s[is].user_info<fastjet::MyUserInfo>().q();
          int matching_s = jets_s[is].user_info<fastjet::MyUserInfo>().matching();
          jets_s[is].set_user_info(new fastjet::MyUserInfo(jet_q_s, matching_s, (*p)->pdg_id()));
          break;
        }
      }

      // Match smeared jets with unsmeared jets
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        for (unsigned int is = 0; is < jets_s.size(); is++)
        {
          if (jets_s[is].user_info<fastjet::MyUserInfo>().matching() == matched)
            continue;
          if (jets[i].delta_R(jets_s[is]) < 0.4)
          {
            vector<double> jet_q = jets[i].user_info<fastjet::MyUserInfo>().q();
            vector<double> jet_q_s = jets_s[is].user_info<fastjet::MyUserInfo>().q();
            int parton = jets[i].user_info<fastjet::MyUserInfo>().parton();
            int parton_s = jets_s[is].user_info<fastjet::MyUserInfo>().parton();
            jets[i].set_user_info(new fastjet::MyUserInfo(jet_q, matched, parton));
            jets_s[is].set_user_info(new fastjet::MyUserInfo(jet_q_s, matched, parton_s));
            for (vector<double>::iterator it = jet_q.begin(); it != jet_q.end(); ++it)
            {
              mytxtfile << *it << ", ";
            }
            for (vector<double>::iterator it = jet_q_s.begin(); it != jet_q_s.end(); ++it)
            {
              mytxtfile << *it << ", ";
            }
            mytxtfile << xsecweight << ", " << matched << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << parton << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets[i].constituents().size() << ", " << jets_s[is].perp() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << parton_s << ", " << sdjets_s[is].m() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets_s[is].constituents().size() << ", " << i << ", " << jets[i].delta_R(partons.at(0)) << ", " << jets[i].delta_R(partons.at(1)) << ", " << theEvent->event_scale() << "\n";
            break;
          }
        }
      }

      for (unsigned int i = 0; i < jets.size(); i++)
      {
        if (jets[i].user_info<fastjet::MyUserInfo>().matching() == matched)
          continue;
        vector<double> jet_q = jets[i].user_info<fastjet::MyUserInfo>().q();
        int parton = jets[i].user_info<fastjet::MyUserInfo>().parton();
        jets[i].set_user_info(new fastjet::MyUserInfo(jet_q, missing, parton));

	for (vector<double>::iterator it = jet_q.begin(); it != jet_q.end(); ++it)
        {
          mytxtfile << *it << ", -999, ";
        }
        mytxtfile << xsecweight << ", " << missing << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << parton << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets[i].constituents().size() << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << i << ", " << jets[i].delta_R(partons.at(0)) << ", " << jets[i].delta_R(partons.at(1)) << ", " << theEvent->event_scale() << "\n";
      }
      for (unsigned int is = 0; is < jets_s.size(); is++)
      {
        if (jets_s[is].user_info<fastjet::MyUserInfo>().matching() == matched)
          continue;
        vector<double> jet_q_s = jets_s[is].user_info<fastjet::MyUserInfo>().q();
        int parton_s = jets_s[is].user_info<fastjet::MyUserInfo>().parton();
        jets_s[is].set_user_info(new fastjet::MyUserInfo(jet_q_s, fake, parton_s));

 	mytxtfile << "-999, ";
        for (vector<double>::iterator it = jet_q_s.begin(); it != jet_q_s.end(); ++it)
        {
          mytxtfile << *it << ", ";
        }
        mytxtfile << xsecweight << ", " << fake << ", -999, -999, -999, -999, -999, -999, -999, -999, " << jets_s[is].perp() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << parton_s << ", " << sdjets_s[is].m() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets_s[is].constituents().size() << ", -999, -999, -999, -999\n";
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

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_0810);
}
