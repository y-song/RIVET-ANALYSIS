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
  class ParticleInfo : public PseudoJet::UserInfoBase
  {
  public:
    ParticleInfo(const double charge, const int &pid, const bool &is_direct, const double &tprod, const vector<Rivet::Particle> &parents) : _charge(charge), _pid(pid), _is_direct(is_direct), _tprod(tprod), _parents(parents) {}
    double charge() const { return _charge; }
    int pid() const { return _pid; }
    bool is_direct() const { return _is_direct; }
    double tprod() const { return _tprod; }
    vector<Rivet::Particle> parents() const { return _parents; }

  protected:
    double _charge;
    int _pid;
    bool _is_direct;
    double _tprod;
    vector<Rivet::Particle> _parents;
  };
}

namespace Rivet
{
  class PYTHIA_JETS_RC_FOR_CHARLES : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_RC_FOR_CHARLES);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");
      // open output file
      mytxtfile.open("pythia1208_PTMINpthatPTMAX_rc_truth_tprod.txt");
      // evttxtfile.open("pythia_PTMINpthatPTMAX_omega.txt");
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
        if (abs(p.pid()) == 12 || abs(p.pid()) == 14 || abs(p.pid()) == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        int pid = p.pid();
        bool is_direct = p.isDirect();
        double tprod = p.origin().t() * pow(10, 12); // fm = 10^12 mm
        /* this is equivalent to the following:
        HepMC::GenParticle *theParticle = (HepMC::GenParticle *)p.genParticle();
        double tprod = theParticle->production_vertex()->position()->t(); 
        But it doesn't give the tProd() value that PYTHIA gives,
        since the distribution is different from slide 12 of http://home.thep.lu.se/~torbjorn/talks/cern21accuracy.pdf,
        especially since my code gives 0 for direct particles but tProd() from PYTHIA doesn't */
        Particles parents = p.parents(Cuts::abspid > 10 && Cuts::abspid != 21); // don't save parents if they are partons
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid, is_direct, tprod, parents));
        parts.push_back(pseudojet);
      }

      //===========================================
      // process jets
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::contrib::SoftDrop sd(0, 0.1);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
      // vector<fastjet::PseudoJet> sdjets = sd(jets);
      vector<fastjet::PseudoJet> partons = {};

      for (unsigned int i = 0; i < jets.size(); i++)
      {
        jets[i].set_user_index(0);
      }

      for (unsigned int i = 0; i < jets.size(); i++) // jet loop
      {

        vector<fastjet::PseudoJet> IncParts = sorted_by_pt(jets[i].constituents());
        if (IncParts.size() < 2)
          continue;
        PseudoJet leading_part = IncParts.at(0);
        PseudoJet subleading_part = IncParts.at(1);

        // for charles, skip jets that don't have charged leading or subleading
        if (leading_part.user_info<fastjet::ParticleInfo>().charge() == 0 || subleading_part.user_info<fastjet::ParticleInfo>().charge() == 0)
          continue;

        int parent1 = -999;
        int parent2 = -999;
        Particles parents1 = leading_part.user_info<fastjet::ParticleInfo>().parents();
        Particles parents2 = subleading_part.user_info<fastjet::ParticleInfo>().parents();
        if (parents1.size() > 0)
        {
          parent1 = parents1.at(0).pid();
        }
        if (parents2.size() > 0)
        {
          parent2 = parents2.at(0).pid();
        }

        if (parents1.size() > 1) // this should not happen
        {
          cout << "\nleading track pid: " << IncParts.at(0).user_info<fastjet::ParticleInfo>().pid() << ", parent pids: ";
          for (Particle parent1 : parents1)
          {
            cout << parent1.pid() << ", ";
          }
        }

        mytxtfile << evid << ", " << xsecweight << ", " << jets[i].perp() << ", " << jets[i].constituents().size() << ", " << jets[i].user_index() << ", " << IncParts.at(0).perp() << ", " << IncParts.at(0).user_info<fastjet::ParticleInfo>().charge() << ", " << IncParts.at(0).eta() << ", " << IncParts.at(0).user_info<fastjet::ParticleInfo>().is_direct() << ", " << IncParts.at(0).user_info<fastjet::ParticleInfo>().pid() << ", " << IncParts.at(0).user_info<fastjet::ParticleInfo>().tprod() << ", " << parent1 << ", " << IncParts.at(1).perp() << ", " << IncParts.at(1).user_info<fastjet::ParticleInfo>().charge() << ", " << IncParts.at(1).eta() << ", " << IncParts.at(1).user_info<fastjet::ParticleInfo>().is_direct() << ", " << IncParts.at(1).user_info<fastjet::ParticleInfo>().pid() << ", " << IncParts.at(1).user_info<fastjet::ParticleInfo>().tprod() << ", " << parent2 << "\n";
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      double norm = crossSection() / sumOfWeights();
      std::cout << norm << std::endl;
      mytxtfile.close();
      // evttxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
    // std::ofstream evttxtfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_RC_FOR_CHARLES);
}
