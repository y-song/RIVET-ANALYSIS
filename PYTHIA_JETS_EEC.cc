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
    MyUserInfo(const double &q, const int &matching, const int &parton) : _q(q), _matching(matching), _parton(parton) {}
    double q() const { return _q; }
    int matching() const { return _matching; }
    int parton() const { return _parton; }

  protected:
    double _q;
    int _matching;
    int _parton;
  };

  class ParticleInfo : public PseudoJet::UserInfoBase
  {
  public:
    ParticleInfo(const double charge, const int &pid, const int &baryon) : _charge(charge), _pid(pid), _baryon(baryon) {}
    double charge() const { return _charge; }
    int pid() const { return _pid; }
    int baryon() const { return _baryon; }

  protected:
    double _charge;
    int _pid;
    int _baryon;
  };
}

namespace Rivet
{
  class PYTHIA_JETS_EEC : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_EEC);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.); // 2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      trkfnfile.open("pythia_PTMINpthatPTMAX_trkfn_kp.txt");
      eecfile.open("pythia_PTMINpthatPTMAX_eec_kp.txt");
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      double xsecweight = handler().nominalCrossSection();
      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();
     
      //===========================================
      // process particles
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();

      PseudoJets parts;

      for (const Particle &p : fsParticles)
      {
        if (p.pid() == 12 || p.pid() == 14 || p.pid() == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        int pid = p.pid();
        int baryon = 0;
        if (p.isBaryon() && pid < 0)
        {
          baryon = -1;
        }
        else if (p.isBaryon() && pid > 0)
        {
          baryon = 1;
        }
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid, baryon));
        parts.push_back(pseudojet);
      }

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(600.0) * fastjet::SelectorEtaRange(-2.0, 2.0);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));

      //-------------------------------------------
      // match jets to outgoing partons
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        jets[i].set_user_info(new fastjet::MyUserInfo(-999, -999, -999));
      }

      vector<fastjet::PseudoJet> partons = {};
      for (HepMC::GenEvent::particle_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
      {
        if ((*p)->status() != 23)
          continue;
        PseudoJet parton_pseudojet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        partons.push_back(parton_pseudojet);

        for (unsigned int i = 0; i < jets.size(); i++)
        {
          if (jets[i].delta_R(parton_pseudojet) > 0.4 || jets[i].user_info<fastjet::MyUserInfo>().parton() != -999)
            continue; // skip if delta R is large OR the jet has already been matched
          jets[i].set_user_info(new fastjet::MyUserInfo(-999, -999, (*p)->pdg_id()));
          break;
        }
      }

      //-------------------------------------------
      // loop over all jets in the event
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        double xkplus = 0.;
        double xkminus = 0.;
        double xpiplus = 0.;
        double xpiminus = 0.;
        int nkplus = 0;
        int nkminus = 0;
        int npiplus = 0;
        int npiminus = 0;

        // constituent loop
        for (unsigned int j = 0; j < jets[i].constituents().size(); j++)
        {
          PseudoJet part = jets[i].constituents().at(j);
          if (part.user_info<fastjet::ParticleInfo>().pid() == 321)
          {
            nkplus += 1;
            xkplus += part.perp() / jets[i].perp();
          }
          else if (part.user_info<fastjet::ParticleInfo>().pid() == -321)
          {
            nkminus += 1;
            xkminus += part.perp() / jets[i].perp();
          }
          else if (part.user_info<fastjet::ParticleInfo>().pid() == 211)
          {
            npiplus += 1;
            xpiplus += part.perp() / jets[i].perp();
          }
          else if (part.user_info<fastjet::ParticleInfo>().pid() == -211)
          {
            npiminus += 1;
            xpiminus += part.perp() / jets[i].perp();
          }
        }
        trkfnfile << evid << ", " << i << ", " << xsecweight << ", " << theEvent->event_scale() << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << nkplus << ", " << nkminus << ", " << npiplus << ", " << npiminus << ", " << xkplus << ", " << xkminus << ", " << xpiplus << ", " << xpiminus << ", " << jets[i].user_info<fastjet::MyUserInfo>().parton() << endl;

        // loop over all kaons in the jet
        for (unsigned int k = 0; k < jets[i].constituents().size(); k++)
        {
          PseudoJet kaon_candidate = jets[i].constituents().at(k);
          if (abs(kaon_candidate.user_info<fastjet::ParticleInfo>().pid()) != 321)
          {
            continue;
          }
          // loop over all kaons+pions in the jet
          for (unsigned int p = 0; p < jets[i].constituents().size(); p++)
          {
            if (p == k)
            {
              continue;
            }
            PseudoJet assoc = jets[i].constituents().at(p);
            if (abs(assoc.user_info<fastjet::ParticleInfo>().pid()) != 211 && abs(assoc.user_info<fastjet::ParticleInfo>().pid()) != 321)
            {
              continue;
            }

            int pid_sum = kaon_candidate.user_info<fastjet::ParticleInfo>().pid() + assoc.user_info<fastjet::ParticleInfo>().pid();
            int code = -9;
            double weight = kaon_candidate.e() * assoc.e() / pow(jets[i].perp(), 2.0);
            double dr = kaon_candidate.delta_R(assoc);
            if (pid_sum == 0) // k+k-
            {
              code = 0;
            }
            else if (pid_sum == 642) // k+k+
            {
              code = 1;
            }
            else if (pid_sum == -642) // k-k-
            {
              code = -1;
            }
            else if (pid_sum == 532) // k+pi+
            {
              code = 2;
            }
            else if (pid_sum == -532) // k-pi-
            {
              code = -2;
            }
            else if (pid_sum == 110) // k+pi-
            {
              code = 3;
            }
            else if (pid_sum == -110) // k-pi+
            {
              code = -3;
            }
            eecfile << evid << ", " << i << ", " << xsecweight << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << jets[i].user_info<fastjet::MyUserInfo>().parton() << ", " << pid_sum << ", " << weight << ", " << dr << endl;
          } // track pair loop ends
        }   // kaon loop ends
      }     // jet loop ends
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      trkfnfile.close();
      eecfile.close();
    }

  private:
    std::ofstream trkfnfile;
    std::ofstream eecfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_EEC);
}
