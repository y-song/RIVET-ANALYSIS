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
    ParticleInfo(const double charge, const int &pid, const bool &from_hadron, const vector<Rivet::Particle> &parents) : _charge(charge), _pid(pid), _from_hadron(from_hadron), _parents(parents) {}
    double charge() const { return _charge; }
    int pid() const { return _pid; }
    bool from_hadron() const { return _from_hadron; }
    vector<Rivet::Particle> parents() const { return _parents; }

  protected:
    double _charge;
    int _pid;
    bool _from_hadron;
    vector<Rivet::Particle> _parents;
  };
}

namespace Rivet
{
  class PYTHIA_JETS_RC_TRUTH : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_RC_TRUTH);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 2.0); // !!! special setting for eta check!
      declare(fs, "fs");

      // open output file
      mytxtfile.open("pythia_PTMINpthatPTMAX_eta.txt");
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
        if ( abs(p.pid()) == 12 || abs(p.pid()) == 14 || abs(p.pid()) == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        int pid = p.pid();
        bool from_hadron = p.fromHadron();
        Particles parents = p.parents();
	pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid, from_hadron, parents));
        parts.push_back(pseudojet);
      }

      //===========================================
      // process jets
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-1.6, 1.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::contrib::SoftDrop sd(0, 0.1);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));
      //vector<fastjet::PseudoJet> sdjets = sd(jets);

      for (unsigned int i = 0; i < jets.size(); i++)
      {
        
        //vector<fastjet::PseudoJet> IncParts = sorted_by_pt(jets[i].constituents());
        //vector<fastjet::PseudoJet> ChargedParts;

        // loop over all constituents to create a vector of charged constituents
        /*for (unsigned int j = 0; j < IncParts.size(); j++)
        {
          PseudoJet part = IncParts.at(j);
          double charge = part.user_info<fastjet::ParticleInfo>().charge();
          if (charge != 0)
          {
            ChargedParts.push_back(part);
          }
        }

        if (ChargedParts.size() < 2)
          continue;
            
        double dr = ChargedParts.at(0).delta_R(ChargedParts.at(1));
	double m2 = pow(ChargedParts.at(0).e()+ChargedParts.at(1).e(), 2.0) - pow(ChargedParts.at(0).px()+ChargedParts.at(1).px(), 2.0) - pow(ChargedParts.at(0).py()+ChargedParts.at(1).py(), 2.0) - pow(ChargedParts.at(0).pz()+ChargedParts.at(1).pz(), 2.0);
        double m = sqrt(m2);
	//cout << m << endl;

	int parent1 = -999;
	int parent2 = -999;
	Particles parents1 = ChargedParts.at(0).user_info<fastjet::ParticleInfo>().parents();
	Particles parents2 = ChargedParts.at(1).user_info<fastjet::ParticleInfo>().parents();
	if (ChargedParts.at(0).user_info<fastjet::ParticleInfo>().from_hadron())
	{
	    if (parents1.size() != 1) // this only happens rarely
	    {
		cout << "\nleading track pid: " << ChargedParts.at(0).user_info<fastjet::ParticleInfo>().pid() << ", parent pids";
       		for (Particle parent1 : parents1)
		{
	    	    cout << parent1.pid() << ", ";
		}
		continue;
	    }
	    parent1 = parents1.at(0).pid();
	}
	if (ChargedParts.at(1).user_info<fastjet::ParticleInfo>().from_hadron())
	{
	    if (parents2.size() != 1)
		continue;
	    parent2 = parents2.at(0).pid();
	}*/
	mytxtfile << evid << ", " << xsecweight << ", " << jets[i].perp() << ", " << jets[i].eta() << "\n";
	//mytxtfile << evid << ", " << xsecweight << ", " << jets[i].perp() << ", " << jets[i].constituents().size() << ", " << ChargedParts.size() << ", " << ChargedParts.at(0).perp() << ", " << ChargedParts.at(0).user_info<fastjet::ParticleInfo>().charge() << ", " << ChargedParts.at(0).eta() << ", " << ChargedParts.at(0).phi() << ", " << ChargedParts.at(0).user_info<fastjet::ParticleInfo>().pid() << ", " << parent1 << ", " << ChargedParts.at(1).perp() << ", " << ChargedParts.at(1).user_info<fastjet::ParticleInfo>().charge() << ", " << ChargedParts.at(1).eta() << ", " << ChargedParts.at(1).phi() << ", " << ChargedParts.at(1).user_info<fastjet::ParticleInfo>().pid() << ", " << parent2 << ", " << dr << ", " << m << "\n";
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
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_RC_TRUTH);
}
