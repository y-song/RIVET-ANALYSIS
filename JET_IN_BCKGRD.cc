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
    ParticleInfo(const double &charge, const int &j) : _charge(charge), _j(j) {}
    double charge() const { return _charge; }
    int j() const { return _j; }

  protected:
    double _charge;
    int _j; // 1 if particle is from jet, 0 if it is from background
  };

  class DihadronInfo : public PseudoJet::UserInfoBase
  {
  public:
    DihadronInfo(const double &b_ch, const double &b, const double &pt1_ch, const double &pt2_ch, const double &y1_ch, const double &y2_ch, const double &phi1_ch, const double &phi2_ch, const int &j1_ch, const int &j2_ch, const double &pt1, const double &pt2, const double &y1, const double &y2, const double &phi1, const double &phi2, const int &j1, const int &j2) : _b_ch(b_ch), _b(b), _pt1_ch(pt1_ch), _pt2_ch(pt2_ch), _y1_ch(y1_ch), _y2_ch(y2_ch), _phi1_ch(phi1_ch), _phi2_ch(phi2_ch), _j1_ch(j1_ch), _j2_ch(j2_ch), _pt1(pt1), _pt2(pt2), _y1(y1), _y2(y2), _phi1(phi1), _phi2(phi2), _j1(j1), _j2(j2) {}
    double b_ch() const { return _b_ch; }
    double b() const { return _b; }
    double pt1_ch() const { return _pt1_ch; }
    double pt2_ch() const { return _pt2_ch; }
    double y1_ch() const { return _y1_ch; }
    double y2_ch() const { return _y2_ch; }
    double phi1_ch() const { return _phi1_ch; }
    double phi2_ch() const { return _phi2_ch; }
    int j1_ch() const { return _j1_ch; } // 1 if particle is from jet, 0 if it is from background
    int j2_ch() const { return _j2_ch; }
    double pt1() const { return _pt1; }
    double pt2() const { return _pt2; }
    double y1() const { return _y1; }
    double y2() const { return _y2; }
    double phi1() const { return _phi1; }
    double phi2() const { return _phi2; }
    int j1() const { return _j1; }
    int j2() const { return _j2; }

  protected:
    double _b_ch;
    double _b;
    double _pt1_ch;
    double _pt2_ch;
    double _y1_ch;
    double _y2_ch;
    double _phi1_ch;
    double _phi2_ch;
    int _j1_ch;
    int _j2_ch;
    double _pt1;
    double _pt2;
    double _y1;
    double _y2;
    double _phi1;
    double _phi2;
    int _j1;
    int _j2;
  };
}

namespace Rivet
{
  class JET_IN_BCKGRD : public Analysis
  {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(JET_IN_BCKGRD);

    void init()
    {
      const FinalState fs(Cuts::pt > 0.2 && Cuts::pt < 30 && Cuts::abseta < 1.0);
      declare(fs, "fs");

      // open output file
      mytxtfile.open("jet_in_bckgrd_test.txt");
    }

    // Function to sample from the distribution x * exp(-6 * x)
    double random_mt(std::default_random_engine &generator)
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

    double random_rhoa(std::default_random_engine &generator)
    {

      double mean = 19.;
      double stddev = 4.5;

      // Create a normal distribution object with the desired mean and standard deviation
      std::normal_distribution<double> gaussian(mean, stddev);
      double rhoa = gaussian(generator);

      return rhoa;
    }

    vector<double> get_direction_wrong(double jet_px, double jet_py, double jet_pz)
    {
      
      double random1 = (double)rand() / (double)RAND_MAX;
      double random2 = (double)rand() / (double)RAND_MAX;

      double delta_theta = 0.4 * sqrt(random1);
      double phi_prime = 2 * M_PI * random2;
      double jet_p = sqrt(pow(jet_px, 2.) + pow(jet_py, 2) + pow(jet_pz, 2.));
      double jet_theta = acos(jet_pz/jet_p);
      
      double x_comp = 1/tan(delta_theta) * jet_px / jet_p + cos(phi_prime) * cos(jet_theta);
      double y_comp = 1/tan(delta_theta) * jet_py / jet_p + sin(phi_prime);
      double z_comp = 1/tan(delta_theta) * jet_pz / jet_p + cos(phi_prime) * sin(jet_theta);
      double norm = sqrt(pow(x_comp, 2.) + pow(y_comp, 2.) + pow(z_comp, 2.));
      x_comp /= norm;
      y_comp /= norm;
      z_comp /= norm;

      vector<double> direction = {x_comp, y_comp, z_comp};
      return direction;

    }

    void get_direction(double *rad, double *angle) // thanks Andrew!
    {
      double random1 = (double)rand() / (double)RAND_MAX;
      double random2 = (double)rand() / (double)RAND_MAX;

      *rad = 0.4 * sqrt(random1);
      *angle = 2 * M_PI * random2;
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      HepMC::GenEvent *theEvent = (HepMC::GenEvent *)event.genEvent();
      int evid = theEvent->event_number();
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();
      PseudoJets parts;

      for (const Particle &p : fsParticles)
      {
        if (p.pid() == 12 || p.pid() == 14 || p.pid() == 16)
          continue;
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge, 1));
        parts.push_back(pseudojet);
      }

      // cluster particles into jets
      fastjet::Selector selector = fastjet::SelectorPtMin(15.0) * fastjet::SelectorEtaRange(-0.6, 0.6);
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(selector(cs.inclusive_jets()));

      std::random_device rd;
      std::default_random_engine generator(rd());
      double rhoA_threshold = random_rhoa(generator);

      for (unsigned int i = 0; i < jets.size(); i++)
      {

        vector<fastjet::PseudoJet> IncPart_unsort = jets[i].constituents();
        vector<fastjet::PseudoJet> ChargedPart_unsort;
        for (unsigned int j = 0; j < IncPart_unsort.size(); j++)
        {
          PseudoJet part = IncPart_unsort.at(j);
          double charge = part.user_info<fastjet::ParticleInfo>().charge();
          if (charge != 0)
          {
            // part.set_user_info(new fastjet::ParticleInfo(charge, 1));
            ChargedPart_unsort.push_back(part);
          }
        }

        double rhoA = 0;

        while (rhoA < rhoA_threshold)
        {

          std::random_device rd;
          std::default_random_engine generator(rd());
          double rad, angle;
          int charge;

          // set kinematics
          get_direction(&rad, &angle);
          double mt = random_mt(generator) + 0.2;
          double pt = sqrt(pow(mt, 2) - 0.13957 * 0.13957);
          double eta = jets[i].eta() + rad * sin(angle);
          double phi = jets[i].phi() + rad * cos(angle);
          double pz = sinh(eta) * pt;
          double py = sin(phi) * pt;
          double px = cos(phi) * pt;
          double E = sqrt(pow(px, 2) + pow(py, 2) + pow(pz, 2) + 0.13957 * 0.13957);

          fastjet::PseudoJet pseudojet(px, py, pz, E);

          // set charge
          std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
          double random = uniformDistribution(generator);
          if (random < 1 / 3)
          {
            charge = -1;
          }
          else if (random > 1 / 3 && random < 2 / 3)
          {
            charge = 0;
          }
          else
          {
            charge = 1;
          }
          pseudojet.set_user_info(new fastjet::ParticleInfo(charge, 0));

          IncPart_unsort.push_back(pseudojet);
          if (charge != 0)
          {
            ChargedPart_unsort.push_back(pseudojet);
          }
          rhoA += pseudojet.pt();

        }

        /* for testing
        for (unsigned int j = 0; j < IncPart_unsort.size(); j++)
        {
          PseudoJet part = IncPart_unsort[j];
          mytxtfile << evid << ", " << part.pt() << ", " << part.eta() << ", " << part.phi() << ", " << part.user_info<fastjet::ParticleInfo>().j() << endl;
        } */

        vector<fastjet::PseudoJet> IncPart = sorted_by_pt(IncPart_unsort);
        vector<fastjet::PseudoJet> ChargedPart = sorted_by_pt(ChargedPart_unsort);

        double b_ch = 0;
        double b = 0;
        double pt1_ch = 0;
        double pt2_ch = 0;
        double y1_ch = 0;
        double y2_ch = 0;
        double phi1_ch = 0;
        double phi2_ch = 0;
        int j1_ch = -9;
        int j2_ch = -9;
        double pt1 = 0;
        double pt2 = 0;
        double y1 = 0;
        double y2 = 0;
        double phi1 = 0;
        double phi2 = 0;
        int j1 = -9;
        int j2 = -9;

        if (IncPart.size() < 2)
        {
          b = -9;
          pt1 = -9;
          pt2 = -9;
          y1 = -9;
          y2 = -9;
          phi1 = -9;
          phi2 = -9;
        }
        else
        {
          b = (IncPart.at(0)).user_info<fastjet::ParticleInfo>().charge() * (IncPart.at(1)).user_info<fastjet::ParticleInfo>().charge();
          pt1 = (IncPart.at(0)).pt();
          pt2 = (IncPart.at(1)).pt();
          y1 = (IncPart.at(0)).rap();
          y2 = (IncPart.at(1)).rap();
          phi1 = (IncPart.at(0)).phi();
          phi2 = (IncPart.at(1)).phi();
          j1 = (IncPart.at(0)).user_info<fastjet::ParticleInfo>().j();
          j2 = (IncPart.at(1)).user_info<fastjet::ParticleInfo>().j();
        }

        if (ChargedPart.size() < 2)
        {
          b_ch = -9;
          pt1_ch = -9;
          pt2_ch = -9;
          y1_ch = -9;
          y2_ch = -9;
          phi1_ch = -9;
          phi2_ch = -9;
        }
        else
        {
          b_ch = (ChargedPart.at(0)).user_info<fastjet::ParticleInfo>().charge() * (ChargedPart.at(1)).user_info<fastjet::ParticleInfo>().charge();
          pt1_ch = (ChargedPart.at(0)).pt();
          pt2_ch = (ChargedPart.at(1)).pt();
          y1_ch = (ChargedPart.at(0)).rap();
          y2_ch = (ChargedPart.at(1)).rap();
          phi1_ch = (ChargedPart.at(0)).phi();
          phi2_ch = (ChargedPart.at(1)).phi();
          j1_ch = (ChargedPart.at(0)).user_info<fastjet::ParticleInfo>().j();
          j2_ch = (ChargedPart.at(1)).user_info<fastjet::ParticleInfo>().j();
        }
        jets[i].set_user_info(new fastjet::DihadronInfo(b_ch, b, pt1_ch, pt2_ch, y1_ch, y2_ch, phi1_ch, phi2_ch, j1_ch, j2_ch, pt1, pt2, y1, y2, phi1, phi2, j1, j2));

	mytxtfile << evid << ", " << rhoA_threshold << ", " << i << ", " << jets[i].user_info<fastjet::DihadronInfo>().b_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().b() << ", " << jets[i].user_info<fastjet::DihadronInfo>().pt1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().pt2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().j1_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().j2_ch() << ", " << jets[i].user_info<fastjet::DihadronInfo>().pt1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().pt2() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().y2() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().phi2() << ", " << jets[i].user_info<fastjet::DihadronInfo>().j1() << ", " << jets[i].user_info<fastjet::DihadronInfo>().j2() << ", " << jets.size() << ", " << jets[i].perp() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << jets[i].constituents().size() << ", " << IncPart.size() << endl;
	
	IncPart.clear();
        ChargedPart.clear();
        IncPart_unsort.clear();
        ChargedPart_unsort.clear();

      } // jet loop ends
   } // event loop ends

    /// Normalise histograms etc., after the run
    void finalize()
    {
      mytxtfile.close();
    }

  private:
    std::ofstream mytxtfile;
  };

  RIVET_DECLARE_PLUGIN(JET_IN_BCKGRD);
}
