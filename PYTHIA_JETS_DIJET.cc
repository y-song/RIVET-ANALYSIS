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
  class PYTHIA_JETS_DIJET : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_DIJET);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      eecfile.open("pythia_PTMINpthatPTMAX_eec.txt");
    }

    float angle(fastjet::PseudoJet &a, fastjet::PseudoJet &b)
    {
	double ax = a.px();
    	double ay = a.py();
    	double az = a.pz();
    	double bx = b.px();
    	double by = b.py();
    	double bz = b.pz();

    	double dot = ax * bx + ay * by + az * bz;
    	double mag_a = sqrt(ax * ax + ay * ay + az * az);
    	double mag_b = sqrt(bx * bx + by * by + bz * bz);
    	double cos_theta = dot / (mag_a * mag_b);
    	return acos(cos_theta);
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
        PseudoJet pseudojet( p.px(), p.py(), p.pz(), p.E() ); // assume masssive particles
        double charge = p.charge();
        int pid = p.pid();
        int baryon = 0;
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, pid, baryon));
        parts.push_back(pseudojet);
      }

      //===========================================
      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(10.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));

      //-------------------------------------------
      // match jets to outgoing partons
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        jets[i].set_user_info(new fastjet::MyUserInfo(-999, -999, -999));
      }

      /*vector<fastjet::PseudoJet> partons = {};
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
      }*/

      //-------------------------------------------
      // loop over leading dijets in the event
      if (jets.size() >= 2 && jets[0].perp() >= 15)
      {
    
      for (unsigned int i = 0; i < 2; i++)
      {
       
        // loop over all particles in jet i
        for (unsigned int k = 0; k < jets[i].constituents().size(); k++)
        {
          PseudoJet particle1 = jets[i].constituents().at(k);
          if (particle1.user_info<fastjet::ParticleInfo>().charge() == 0)
          {
	    continue;
	  }

          for (unsigned int j = i; j < 2; j++)
		  
	  {
	    // loop over all particles in jet j
            for (unsigned int p = 0; p < jets[j].constituents().size(); p++)
            {
              if (p == k && i==j)
              {
                continue;
              }
              PseudoJet particle2 = jets[j].constituents().at(p); 
	      if (particle2.user_info<fastjet::ParticleInfo>().charge() == 0)
	      {
	        continue;
	      }
  
	      double weight = particle1.e() * particle2.e() / pow((jets[i].perp()+jets[j].perp()), 2.0);
    	      //double dr = particle1.delta_R(particle2); // sqrt(delta y^2 + delta phi^2)
	      double dr = angle(particle1, particle2);
  	      double jet_dr = jets[i].delta_R(jets[j]);
  	      double jet_dphi = jets[i].delta_phi_to(jets[j]);

  	      if (dr > 3.14159265 || dr == 0)
	      {
		cout << dr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	      }

	      eecfile << evid << ", " << i << ", " << xsecweight << ", " << jets[i].perp() << ", " << jets[j].perp() << ", " << jet_dr << ", " << jet_dphi << ", " << weight << ", " << dr << endl;
            } // track pair loop ends
          }   // jet j loop ends
	}   // particle1 loop ends
      }     // jet i loop ends
    }       // close if statement
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      eecfile.close();
    }

  private:
    std::ofstream eecfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_DIJET);
}
