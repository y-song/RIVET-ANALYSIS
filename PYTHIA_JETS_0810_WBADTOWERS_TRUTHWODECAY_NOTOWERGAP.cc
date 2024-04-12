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
  class PYTHIA_JETS_0810_WBADTOWERS_TRUTHWODECAY_NOTOWERGAP : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_0810_WBADTOWERS_TRUTHWODECAY_NOTOWERGAP);

    void init()
    {
      const FinalState fs(Cuts::open());
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
			mytxtfile.open("pythia_5pthat9_0923.txt");

      // make JP grids
      cout << "eta bins: " << endl;
      for (int i = 0; i <= Nbounds_JP_eta; ++i)
      {
        double etax = -1 * etaMax + i * (double)2 * etaMax / Nbounds_JP_eta;
        etabins_JP.push_back(etax);
        cout << etax << endl;
      }
      cout << "phi bins: " << endl;
      for (int i = 0; i <= Nbounds_JP_phi; ++i)
      {
        double phix = i * (double)phiMax / Nbounds_JP_phi;
        phibins_JP.push_back(phix);
        cout << phix << endl;
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
      int reject = 0;
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

      //-------------------------------------------
      // fill pythia-level particles (undo decay)
      Particles pythiaParticles;
      Particles all_parents;
      vector<int> all_children_index; // index in fsParticles vector

      pythiaParticles.clear();
      all_parents.clear();
      all_children_index.clear();

      for (unsigned int i = 0; i < fsParticles.size(); i++)
      {
        Particle p = fsParticles.at(i);
        vector<Particle> parents = p.parents();
        pythiaParticles.push_back(p); // copy fsParticles to pythiaParticles
        Particle parent = parents[0];
        // check if the particle is from decay
        if (p.fromHadron() && (parent.abspid() == 111 || parent.abspid() == 211 || parent.abspid() == 221 || parent.abspid() == 321 || parent.abspid() == 310 || parent.abspid() == 130 || parent.abspid() == 3122 || parent.abspid() == 3212 || parent.abspid() == 3112 || parent.pid() == 3222 || parent.abspid() == 3312 || parent.abspid() == 3322 || parent.abspid() == 3334)) // from a hadron weak decay
        {
          all_children_index.push_back(i);
          if (parents.size() > 1)
          {
            cout << "number of parents: " << parents.size() << endl;
          }
          // check if the parent is already in the vector
          bool recorded = false;
          for (const Particle &it : all_parents)
          {
            if (parent.isSame(it))
            {
              recorded = true;
              break;
            }
          }
          if (recorded)
            continue;
          all_parents.push_back(parent);
        }
      }
      //  remove all_children from pythiaParticles by index
      if (all_children_index.size() != 0)
      {
        for (int i = all_children_index.size() - 1; i > -1; i--)
        {
          if (pythiaParticles.size() > 0)
          {
            pythiaParticles.erase(pythiaParticles.begin() + all_children_index[i]);
          }
          else
          {
            cout << all_children_index[i] << endl;
          }
        }
      }
      //  include all_parents in pythiaParticles
      for (const Particle &p : all_parents)
      {
        pythiaParticles.push_back(p);
      }

      for (const Particle &p : pythiaParticles)
      {
        if (p.abseta() > 1.0 || p.pt() < 0.2 || p.pt() > 30)
          continue;
        if (p.abspid() == 12 || p.abspid() == 14 || p.abspid() == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        int pid = p.pid();
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid));
        parts.push_back(pseudojet);
      }

      //-------------------------------------------
      // fill geant-level particles (keep decay products)
      for (const Particle &p : fsParticles)
      {
        if (p.abspid() == 12 || p.abspid() == 14 || p.abspid() == 16)
          continue;
        double charge = p.charge();
        int pid = p.pid();
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
          if (pj_smeared.Et() > 30.)
          {
            reject = 3;
            continue;
          }
          if (pj_smeared.Et() < 0.2 || abs(pj_smeared.eta()) > 1.0)
            continue;
          bool bad_tower = false;
          double phi_smeared;
          if (pj_smeared.phi() > pi)
          {
            phi_smeared = pj_smeared.phi() - 2.0 * pi;
          }
          else if (pj_smeared.phi() < -1.0 * pi)
          {
            phi_smeared = pj_smeared.phi() + 2.0 * pi;
          }
          else
          {
            phi_smeared = pj_smeared.phi();
          }
          for (int j = 0; j < 235; j++)
          {
            if ((abs(pj_smeared.eta()) - bad_eta_list[j]) < 0.025 && abs(phi_smeared - bad_phi_list[j]) < 0.025)
            // neutral particle eta and phi fall in bad eta and phi lists
            {
              bad_tower = true;
              break;
            }
          }
          if (bad_tower)
          {
            continue;
          }
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
          if (pj_smeared.perp() > 30.)
          {
            reject = 3;
            continue;
          }
          if (pj_smeared.perp() < 0.2 || abs(pj_smeared.eta()) > 1.0)
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
          //cout << x << endl;
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
            //cout << "jp2_eta, jp2_phi: " << x << ", " << y << endl;
            break;
          }
        }
      }

      //===========================================
      // jet selectors
      // THIS IS NOT THE DEFAULT!
      fastjet::Selector selector = fastjet::SelectorPtMin(5.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::Selector selector_s = fastjet::SelectorPtMin(0.2) * fastjet::SelectorEtaRange(-0.6, 0.6);

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
      for (unsigned int i = 0; i < jets_s.size(); i++)
      {
        vector<double> q{};
        for (const double &kappa : k)
        {
          double numerator = 0;
          for (unsigned int j = 0; j < jets_s[i].constituents().size(); j++)
          {
            PseudoJet part = jets_s[i].constituents().at(j);
            double charge = part.user_info<fastjet::ParticleInfo>().charge();
            numerator += pow(part.perp(), kappa) * charge;
          }
          q.push_back(numerator / pow(jets_s[i].perp(), kappa));
        }
        jets_s[i].set_user_info(new fastjet::MyUserInfo(q, -1, -999));
      }

      //-------------------------------------------
      // Match jets to partons
      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();

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

      //-------------------------------------------
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
            mytxtfile << xsecweight << ", " << evid << ", " << matched << ", " << reject << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << parton << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << jets_s[is].perp() << ", " << jets_s[is].eta() << ", " << jets_s[is].phi() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << parton_s << ", " << sdjets_s[is].m() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << i << ", " << jets[i].delta_R(partons.at(0)) << ", " << jets[i].delta_R(partons.at(1)) << ", " << theEvent->event_scale() << ", " << isJP2 << ", " << JP2_et << ", " << JP2_eta << ", " << JP2_phi << "\n";
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
        mytxtfile << xsecweight << ", " << evid << ", " << missing << ", " << reject << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << parton << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << i << ", " << jets[i].delta_R(partons.at(0)) << ", " << jets[i].delta_R(partons.at(1)) << ", " << theEvent->event_scale() << ", " << isJP2 << ", " << JP2_et << ", " << JP2_eta << ", " << JP2_phi << "\n";
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
        mytxtfile << xsecweight << ", " << evid << ", " << fake << ", -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, " << jets_s[is].perp() << ", " << jets_s[is].eta() << ", " << jets_s[is].phi() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << parton_s << ", " << sdjets_s[is].m() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", -999, -999, -999, -999" << ", " << isJP2 << ", " << JP2_et << ", " << JP2_eta << ", " << JP2_phi << "\n";
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      double norm = crossSection() / sumOfWeights();
      std::cout << norm << std::endl;
      mytxtfile.close();
    }

    ///@}

  private:
    double eff_array[45][80];
    double bad_eta_list[235] = {0.725, 0.775, 0.325, 0.675, 0.02675, 0.125, 0.325, 0.375, 0.775, 0.02675, 0.125, 0.325, 0.375, 0.675, 0.425, 0.525, 0.225, 0.675, 0.075, 0.775, 0.075, 0.125, 0.725, 0.475, 0.775, 0.825, 0.475, 0.675, 0.875, 0.925, 0.375, 0.525, 0.675, 0.875, 0.575, 0.675, 0.175, 0.475, 0.525, 0.575, 0.625, 0.325, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.075, 0.125, 0.175, 0.225, 0.475, 0.525, 0.575, 0.625, 0.875, 0.075, 0.125, 0.175, 0.225, 0.325, 0.475, 0.525, 0.575, 0.625, 0.875, 0.675, 0.775, 0.125, 0.875, 0.967, 0.175, 0.967, 0.675, 0.725, 0.825, 0.675, 0.325, 0.425, 0.225, 0.325, 0.425, 0.02675, 0.075, 0.02675, 0.275, 0.525, 0.625, 0.02675, 0.875, 0.925, 0.967, 0.02675, 0.225, 0.375, 0.875, 0.925, 0.967, 0.02675, 0.075, 0.125, 0.175, 0.225, 0.675, 0.875, 0.925, 0.02675, 0.075, 0.125, 0.175, 0.225, 0.875, 0.925, 0.967, 0.02675, 0.725, 0.625, 0.425, 0.675, 0.725, 0.775, 0.125, 0.425, 0.375, 0.475, 0.725, 0.425, 0.875, 0.375, 0.725, 0.775, 0.875, 0.625, 0.725, 0.425, 0.475, 0.675, 0.125, 0.275, 0.325, 0.875, 0.925, 0.225, 0.625, 0.175, 0.725, 0.675, 0.875, 0.625, 0.675, 0.875, 0.375, 0.125, 0.425, 0.725, 0.275, 0.625, -0.475, -0.775, -0.967, -0.275, -0.967, -0.475, -0.525, -0.675, -0.625, -0.475, -0.725, -0.075, -0.475, -0.275, -0.875, -0.525, -0.575, -0.325, -0.02675, -0.475, -0.02675, -0.675, -0.675, -0.725, -0.775, -0.425, -0.425, -0.225, -0.575, -0.425, -0.925, -0.967, -0.525, -0.625, -0.925, -0.625, -0.925, -0.075, -0.925, -0.02675, -0.375, -0.275, -0.275, -0.325, -0.675, -0.875, -0.925, -0.967, -0.675, -0.875, -0.967, -0.225, -0.575, -0.175, -0.775, -0.875, -0.925, -0.225, -0.02675, -0.875, -0.925, -0.225, -0.425, 0.02675, 0.675, 0.967, 0.225, 0.325, 0.875, 0.925, 0.475, 0.375, -0.275, -0.375, -0.02675, -0.875, -0.425};
    double bad_phi_list[235] = {1.23186, 1.07198, 1.02242, 1.02242, 0.862539, 0.812977, 0.603537, 0.603537, 0.603537, 0.54838, 0.54838, 0.54838, 0.54838, 0.54838, 0.498817, 0.234221, -0.0247814, -0.0799383, -0.129501, -0.129501, -0.184658, -0.184658, -0.234221, -0.289378, -0.289378, -0.289378, -0.338941, -0.338941, -0.338941, -0.338941, -0.394098, -0.394098, -0.394098, -0.394098, -0.44366, -0.44366, -0.6531, -0.75782, -0.75782, -0.75782, -0.75782, -0.812977, -0.812977, -0.812977, -0.812977, -0.812977, -0.812977, -0.812977, -0.862539, -0.862539, -0.862539, -0.862539, -0.862539, -0.862539, -0.862539, -0.862539, -0.862539, -0.917696, -0.917696, -0.917696, -0.917696, -0.917696, -0.917696, -0.917696, -0.917696, -0.917696, -0.917696, -0.967259, -0.967259, -1.02242, -1.02242, -1.02242, -1.07198, -1.12714, -1.1767, -1.1767, -1.1767, -1.28142, -1.38614, -1.38614, -1.4413, -1.4413, -1.4413, -1.54601, -1.54601, -1.59558, -1.65073, -1.65073, -1.65073, -1.80502, -1.80502, -1.80502, -1.80502, -1.86017, -1.86017, -1.86017, -1.86017, -1.86017, -1.86017, -1.90974, -1.90974, -1.90974, -1.90974, -1.90974, -1.90974, -1.90974, -1.90974, -1.96489, -1.96489, -1.96489, -1.96489, -1.96489, -1.96489, -1.96489, -1.96489, -2.01446, -2.06961, -2.11918, -2.2239, -2.2239, -2.2239, -2.27905, -2.32862, -2.32862, -2.38377, -2.38377, -2.43334, -2.48849, -2.69793, -2.80265, -2.80265, -2.80265, -2.85221, -2.90737, -3.01209, -3.06165, 3.11681, 3.11681, 2.95693, 2.95693, 2.95693, 2.69793, 2.69793, 2.38377, 2.27905, 2.2239, 2.2239, 2.17433, 2.17433, 2.11918, 2.11918, 2.11918, 2.06961, 1.90974, 1.90974, 1.80502, 1.54601, 1.33658, 1.86017, 1.86017, 1.90974, 1.96489, 1.96489, 2.32862, 2.32862, 2.43334, 2.48849, 2.7475, 2.95693, -2.95693, -2.95693, -2.85221, -2.85221, -2.69793, -2.69793, -2.38377, -2.27905, -2.11918, -1.90974, -1.65073, -1.59558, -1.59558, -1.59558, -1.54601, -1.33658, -1.28142, -1.28142, -1.12714, -1.12714, -1.12714, -1.07198, -1.07198, -1.02242, -0.967259, -0.967259, -0.708257, -0.708257, -0.6531, -0.44366, -0.394097, -0.234221, -0.234221, -0.234221, -0.234221, -0.234221, -0.234221, -0.129501, -0.129501, -0.129501, 0.0799385, 0.603537, 0.6531, 0.6531, 0.6531, 0.917697, 0.967259, 1.07198, 1.49086, 1.49086, 1.54602, 1.75545, 0.338941, 0.0247814, -0.708257, -2.06961, -2.11918, -2.17433, -2.43334, 3.11681, 2.27905, 1.96489, -1.80502, -0.967259, 0.289378, 0.498818};
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

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_0810_WBADTOWERS_TRUTHWODECAY_NOTOWERGAP);
}