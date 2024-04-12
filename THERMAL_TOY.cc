#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
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
  class THERMAL_TOY : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(THERMAL_TOY);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      mytxtfile.open("thermal_toy_50pthat_mult125.txt");
    }

    double sampleFromDistribution(std::default_random_engine &generator)
    {
      // Define the distribution parameters
      double lambda = 6.0;

      // Generate a random number between 0 and 1
      std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
      double u = uniformDistribution(generator);

      // Inverse transform sampling
      double x = -std::log(1.0 - u) / lambda;

      return x;
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {
      std::random_device rd;
      std::default_random_engine generator(rd());

      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();
      Particles fsParticles = sortByPt(applyProjection<FinalState>(event, "fs").particles()); // usage from https://gitlab.com/hepcedar/rivet/-/blob/abuckley/analyses/pluginMC/MC_HFJETS.cc#L54
      // int mult = fsParticles.size();

      // specify event multiplicity
      normal_distribution<double> random_mult(125, 10); // gaussian with a mean of 125 and width of 10
      int mult = round(random_mult(rd));
      if (mult <= 0)
      {
        mult = 1;
      }
      if (fsParticles.size() == 0)
      {
        mult = 0;
      }
      // create particles
      vector<fastjet::PseudoJet> parts;
      //exponential_distribution<double> random_pt(6);
      uniform_real_distribution<double> random_eta(-1.0, 1.0);
      uniform_real_distribution<double> random_phi(0, 2 * 3.14159);
      for (int i = 0; i < mult - 1; i++)
      {
        fastjet::PseudoJet part = fastjet::PseudoJet();

        if (i == 0)
        {
          Particle leading = fsParticles.at(0);
          part.reset_PtYPhiM(leading.pt(), leading.eta(), leading.phi(), leading.mass()); // use the kinematics of the leading particle of the event
        }
        else
        {
          part.reset_PtYPhiM(sampleFromDistribution(generator) + 0.2, random_eta(rd), random_phi(rd), 0.); // require other particles to have pT > 0.2 GeV
        }

        //cout << part.pt() << ", " << part.eta() << endl;
        parts.push_back(part);
      }

      // cluster particles into jets
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));

      fastjet::contrib::SoftDrop sd(0, 0.1);
      vector<fastjet::PseudoJet> sdjets = sd(jets);

      for (unsigned int i = 0; i < jets.size(); i++)
      {
        mytxtfile << evid << ", " << i << ", " << mult << ", " << jets.size() << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << sdjets[i].m() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R() << ", " << sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry() << ", " << sdjets[i].constituents().size() << endl;
      }
    } // event loop ends

    /// Normalise histograms etc., after the run
    void finalize()
    {
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(THERMAL_TOY);
}
