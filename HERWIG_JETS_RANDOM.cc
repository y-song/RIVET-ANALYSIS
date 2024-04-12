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
  class HERWIG_JETS_RANDOM : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERWIG_JETS_RANDOM);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      mytxtfile.open("herwig_PTMINpthatPTMAX_random_kt.txt");
    }
    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      double xsecweight = handler().nominalCrossSection();

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
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid));
        parts.push_back(pseudojet);
      }

      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.6, 0.6);

      fastjet::JetDefinition jet_def(fastjet::kt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));

      fastjet::contrib::SoftDrop sd(0, 0.1);
      vector<fastjet::PseudoJet> sdjets = sd(jets);

      fastjet::JetDefinition newjetdef(fastjet::antikt_algorithm, 2.0);
      
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        mytxtfile << evid << ", " << i << ", " << xsecweight << ", " << theEvent->event_scale() << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets[i].constituents().size() << ", ";

        // shuffle jet constituents' pT to break angular ordering
        vector<double> pt_vec;
	vector<PseudoJet> newparts;

        for (unsigned int j = 0; j < jets[i].constituents().size(); j++)  
        {
          PseudoJet part = jets[i].constituents().at(j);
          pt_vec.push_back(part.perp());
        }

        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(pt_vec.begin(), pt_vec.end(), g);
        for (unsigned int j = 0; j < pt_vec.size(); j++)
        {
          PseudoJet part = jets[i].constituents().at(j);
          double ratio = pt_vec.at(j)/part.perp(); // new pT over old pT
	  double newpx = ratio * part.px();
	  double newpy = ratio * part.py();
	  double newe = sqrt( pow(newpx,2.0) + pow(newpy,2.0) + pow(part.pz(),2.0) + pow(part.m(),2.0) );
	  PseudoJet newpart = PseudoJet( newpx, newpy, part.pz(), newe );
          newparts.push_back(newpart);
	}

	// recluster shuffled particles into a jet
	fastjet::ClusterSequence newcs(newparts, newjetdef);
      	vector<PseudoJet> newjets = newcs.inclusive_jets();

       	fastjet::PseudoJet newjet = newjets.at(0);
        fastjet::PseudoJet newsdjet = sd(newjet);
	
	mytxtfile << newjet.perp() << ", " << newjet.eta() << ", " << newjet.phi() << ", " << newjet.m() << ", " << newjet.constituents().size() << ", " << newsdjet.m() << ", " << newsdjet.structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << newsdjet.structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << newsdjet.constituents().size() << endl;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(HERWIG_JETS_RANDOM);
}
