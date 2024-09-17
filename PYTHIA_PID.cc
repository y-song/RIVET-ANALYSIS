// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/UnstableParticles.hh"
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
  class PYTHIA_PID : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_PID);

    void init()
    {
      const UnstableParticles ufs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0 && Cuts::pid==223);
      declare(ufs, "ufs");
     // open output file
      mytxtfile.open("herwig_PTMINpthatPTMAX_omega.txt");
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {
      
      double xsecweight = handler().nominalCrossSection();
      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();

      //===========================================
      // process particles
      Particles ufsParticles = applyProjection<UnstableParticles>(event, "ufs").particles();
      //PseudoJets parts;
      for (const Particle &p : ufsParticles)
      {
        if ( abs(p.pid()) == 223 ) // omega
	{
	    mytxtfile << evid << ", " << xsecweight << ", " << p.perp() << "\n";
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      double norm = crossSection() / sumOfWeights();
      std::cout << norm << std::endl;
      mytxtfile.close();
      //evttxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
    //std::ofstream evttxtfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_PID);
}
