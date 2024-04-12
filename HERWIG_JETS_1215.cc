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
  class HERWIG_JETS_1215 : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERWIG_JETS_1215);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.);//2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      mytxtfile.open("pythia_11pthat15_trkfn_kp.txt");

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
        int baryon = 0;
	if ( p.isBaryon() && pid < 0)
	{
	  baryon = -1;
	}
      	else if ( p.isBaryon() && pid > 0 )
	{
	  baryon = 1;
        }	  
	pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid, baryon));
        parts.push_back(pseudojet);
      }

      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-2.0, 2.0);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));

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
	    	if ( part.user_info<fastjet::ParticleInfo>().pid() == 321 )
	    	{	
			nkplus += 1;
			xkplus += part.perp()/jets[i].perp();
   		}
		else if ( part.user_info<fastjet::ParticleInfo>().pid() == -321 )
	    	{	
			nkminus += 1;
			xkminus += part.perp()/jets[i].perp();
   		}
		else if ( part.user_info<fastjet::ParticleInfo>().pid() == 211 )
	    	{	
			npiplus += 1;
			xpiplus += part.perp()/jets[i].perp();
   		}
		else if ( part.user_info<fastjet::ParticleInfo>().pid() == -211 )
	    	{	
			npiminus += 1;
			xpiminus += part.perp()/jets[i].perp();
   		}
	   }
	   mytxtfile << evid << ", " << i << ", " << xsecweight << ", " << theEvent->event_scale() << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << nkplus << ", " << nkminus << ", " << npiplus << ", " << npiminus << ", " << xkplus << ", " << xkminus << ", " << xpiplus << ", " << xpiminus << endl;		
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

  RIVET_DECLARE_PLUGIN(HERWIG_JETS_1215);
}
