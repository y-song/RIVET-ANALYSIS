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
    vector <double> q() const { return _q; }
    int matching() const { return _matching; }
    int parton() const { return _parton; }

  protected:
    vector <double> _q;
    int _matching;
    int _parton;
  };

  class ParticleInfo : public PseudoJet::UserInfoBase
  {
  public:
    ParticleInfo(const double charge) : _charge(charge) {}
    double charge() const { return _charge; }

  protected:
    double _charge;
  };  
}

namespace Rivet
{
  /// @brief Add a short analysis description here
  class PYTHIA_JETS_CH : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PYTHIA_JETS_CH);

    void init()
    {
      const FinalState fs(Cuts::pt > 0 && Cuts::abseta < 1.0 && Cuts::charge != 0);
      declare(fs, "fs");

      efftxtfile.open("eff.txt");
      if (!efftxtfile.good()) {
         fputs("efficiency file not found\n", stderr); 
	 abort();
      }
      mytxtfile.open("pythia_15_pthat_jets_ch.txt");

      string line;
      for (int i = 0; getline(efftxtfile, line); ++i){
	     istringstream iss (line);
	     string delimiter = ", ";
	     string token;
	     size_t pos = 0;
	     int j = 0;
	     while ((pos = line.find(delimiter)) != string::npos){
		     token = line.substr(0, pos);
		     eff_array[i][j] = stod(token);
		     line.erase(0, pos + delimiter.length());
		     j += 1;
	     }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      int unmatched = -1;
      int matched = 1;
      int missing = 2;
      int fake = 3;
      vector<double> k{0.3, 0.5, 1, 2};
      double xsecweight = handler().nominalCrossSection();

      // Get the final state particles and start clustering jets
      Particles fsParticles = applyProjection<FinalState>(event, "fs").particles();

      // jet selectors
      fastjet::Selector selector = fastjet::SelectorPtMin(10.0)*fastjet::SelectorRapRange(-0.6,0.6);

      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);


      fastjet::contrib::SoftDrop sd(0, 0.1);

      PseudoJets parts, parts_s;

      for (const Particle &p : fsParticles)
      {
        PseudoJet pseudojet(p.px(), p.py(), p.pz(), p.E());
        double charge = p.charge();
        pseudojet.set_user_info(new fastjet::ParticleInfo(charge));
        parts.push_back(pseudojet);

        random_device generator;
        bool neutral = !p.isCharged();
        if (neutral)
        {
          double es, ratio;
          normal_distribution<double> smearNeutral(p.E(), 0.14 * p.E());
          es = smearNeutral(generator);
          if (es < 0.2 || es > 30.)
            continue;
          ratio = es / p.E();
          PseudoJet pj_smeared(p.px() * ratio, p.py() * ratio, p.pz() * ratio, es);
          pj_smeared.reset_momentum(pj_smeared.px(), pj_smeared.py(), pj_smeared.pz(), sqrt(pow(pj_smeared.px(),2)+pow(pj_smeared.py(),2)+pow(pj_smeared.pz(),2)));
	  if (pj_smeared.perp() < 0.2 || pj_smeared.perp() > 30.)
            continue;
          pj_smeared.set_user_info(new fastjet::ParticleInfo(charge));
          parts_s.push_back(pj_smeared);
        }
        else
        {
          uniform_real_distribution<double> dropCharge(0.0, 1.0);
          double keep = dropCharge(generator);
	  double eff;
	  if (p.pt() < 0.2){
		  eff = 0;
	  }
	  else if (p.pt() < 2){
		  eff = eff_array[(int)floor((p.pt()-0.2)/0.04)][(int)floor((p.eta()+1.975)/0.049375)];
	  }
	  else{
		  eff = eff_array[44][(int)floor((p.eta()+1.975)/0.049375)];
	  }
	  if (keep > eff)
	    continue;
	  double pts, ratio;
	  normal_distribution<double> smearCharge(p.pt(), -0.026+0.02*p.pt()+0.003*pow(p.pt(), 2));
	  pts = smearCharge(generator);
	  ratio = pts / p.pt();
	  PseudoJet pj_smeared(p.px() * ratio, p.py() * ratio, p.pz() * ratio, p.E() * ratio);
	  pj_smeared.reset_momentum(pj_smeared.px(), pj_smeared.py(), pj_smeared.pz(), sqrt(pow(pj_smeared.px(),2)+pow(pj_smeared.py(),2)+pow(pj_smeared.pz(),2)+pow(0.1395704,2))); // charged pion mass
	  if (pj_smeared.perp() < 0.2 || pj_smeared.perp() > 30.)	continue;
	  pj_smeared.set_user_info(new fastjet::ParticleInfo(charge));
          parts_s.push_back(pj_smeared);
        }
      }

      fastjet::ClusterSequence cs(parts, jet_def);
      vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      vector<fastjet::PseudoJet> sdjets = sd(jets);
      vector<fastjet::PseudoJet> partons = {};

      fastjet::ClusterSequence cs_s(parts_s, jet_def);
      vector<fastjet::PseudoJet> jets_s = sorted_by_pt(selector(cs_s.inclusive_jets()));
      vector<fastjet::PseudoJet> sdjets_s = sd(jets_s);

      // Calculate jet charge 
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
       		q.push_back( numerator / pow(jets[i].perp(), kappa) );

		numerator = 0;
        	for (unsigned int j = 0; j < jets[i].constituents().size(); j++)
        	{
         	 PseudoJet part = jets[i].constituents().at(j);
         	 if (part.perp() < 2.)	continue;
		 double charge = part.user_info<fastjet::ParticleInfo>().charge();
         	 numerator += pow(part.perp(), kappa) * charge;
       		}
       		q.push_back( numerator / pow(jets[i].perp(), kappa) );
	}
	jets[i].set_user_info(new fastjet::MyUserInfo(q, -1, -999));
        //jets[i].set_user_info(new fastjet::MyUserInfo(q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], -1, -999));
      }

      // Calculate smeared jet charge
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
       		q.push_back( numerator / pow(jets_s[i].perp(), kappa) );

		numerator = 0;
        	for (unsigned int j = 0; j < jets_s[i].constituents().size(); j++)
        	{
         	 PseudoJet part = jets_s[i].constituents().at(j);
         	 if (part.perp() < 2.)	continue;
		 double charge = part.user_info<fastjet::ParticleInfo>().charge();
         	 numerator += pow(part.perp(), kappa) * charge;
       		}
       		q.push_back( numerator / pow(jets_s[i].perp(), kappa) );
	}
	jets_s[i].set_user_info(new fastjet::MyUserInfo(q, -1, -999));
        //jets_s[i].set_user_info(new fastjet::MyUserInfo(q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], -1, -999));
     }
   
      // Match jets to partons
      HepMC::GenEvent *theEvent = (HepMC::GenEvent*) event.genEvent();
      
      for (HepMC::GenEvent::particle_iterator p = theEvent -> particles_begin(); p!= theEvent -> particles_end(); ++p) {
	
	if ( (*p) -> status() != 23 )	continue;
	PseudoJet parton_pseudojet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e()); 
	partons.push_back(parton_pseudojet);
	
	for (unsigned int i = 0; i < jets.size(); i++){
		if (jets[i].delta_R(parton_pseudojet) > 0.4 || jets[i].user_info<fastjet::MyUserInfo>().parton() != -999)	continue;
		vector <double> jet_q = jets[i].user_info<fastjet::MyUserInfo>().q();
		int matching = jets[i].user_info<fastjet::MyUserInfo>().matching();
		jets[i].set_user_info(new fastjet::MyUserInfo(jet_q, matching, (*p) -> pdg_id()));
		break;
	}

	for (unsigned int is = 0; is < jets_s.size(); is++){
		if (jets_s[is].delta_R(parton_pseudojet) > 0.4 || jets_s[is].user_info<fastjet::MyUserInfo>().parton() != -999)	continue;
		vector <double> jet_q_s = jets_s[is].user_info<fastjet::MyUserInfo>().q();
		int matching_s = jets_s[is].user_info<fastjet::MyUserInfo>().matching();
		jets_s[is].set_user_info(new fastjet::MyUserInfo(jet_q_s, matching_s, (*p) -> pdg_id()));
		break;
	}
      }

      // Match smeared jets with unsmeared jets
      for (unsigned int i = 0; i < jets.size(); i++)
      {
        for (unsigned int is = 0; is < jets_s.size(); is++)
        {
          if (jets_s[is].user_info<fastjet::MyUserInfo>().matching() == matched)
            continue;
          if (jets[i].delta_R(jets_s[is]) < 0.4)
          { // found matching true and smeared jets
            vector <double> jet_q = jets[i].user_info<fastjet::MyUserInfo>().q();
            vector <double> jet_q_s = jets_s[is].user_info<fastjet::MyUserInfo>().q();
	    int parton = jets[i].user_info<fastjet::MyUserInfo>().parton();
	    int parton_s = jets_s[is].user_info<fastjet::MyUserInfo>().parton();
            jets[i].set_user_info(new fastjet::MyUserInfo(jet_q, matched, parton));
            jets_s[is].set_user_info(new fastjet::MyUserInfo(jet_q_s, matched, parton_s));
            for (vector<double>::iterator it = jet_q.begin(); it != jet_q.end(); ++it){
		    mytxtfile << *it << ", ";
	    }
	    for (vector<double>::iterator it = jet_q_s.begin(); it != jet_q_s.end(); ++it){
		    mytxtfile << *it << ", ";
	    }
	    double rg, zg, mg, rgs, mgs, zgs;
	    if (sdjets[i] == 0) // did not pass the softdrop requirement
	    {
		    cout << "sdjets[i] == 0" << endl;
		    rg = -2;
		    zg = -2;
		    mg = -2;
	    }
	    else {
		    rg = sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R();
		    zg = sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry();
		    mg = sdjets[i].m();
	    }	    
	    if (sdjets_s[is] == 0) // did not pass the softdrop requirement
	    {
		    rgs = -2;
		    zgs = -2;
		    mgs = -2;
	    }
	    else {
		    rgs = sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R();
		    zgs = sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry();
		    mgs = sdjets_s[is].m();
	    }	    
	    mytxtfile << xsecweight << ", " << matched << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << parton << ", " << mg << ", " << rg << ", " << zg << ", " << jets_s[is].perp() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << jets_s[is].eta() << ", " << jets_s[is].phi() << ", " << parton_s << ", " << mgs << ", " << rgs << ", " << zgs << ", " << i << ", " << jets[i].delta_R(partons.at(0)) << ", " << jets[i].delta_R(partons.at(1)) << ", " << theEvent->event_scale() << "\n";
	    break;
	  }
        }
      }

/*      for (unsigned int i = 0; i < jets.size(); i++)
      {
        if (jets[i].user_info<fastjet::MyUserInfo>().matching() == matched)
          continue;
	// true jets without smeared jets matched
        vector <double> jet_q = jets[i].user_info<fastjet::MyUserInfo>().q();
	int parton = jets[i].user_info<fastjet::MyUserInfo>().parton();
        jets[i].set_user_info(new fastjet::MyUserInfo(jet_q, missing, parton));
        for (vector<double>::iterator it = jet_q.begin(); it != jet_q.end(); ++it){
	    mytxtfile << *it << ", ";
	}
	double rg, zg, mg;
       	if (sdjets[i] == 0) // did not pass the softdrop requirement
	{
	    rg = -2;
	    zg = -2;
	    mg = -2;
	}
	else {
	    rg = sdjets[i].structure_of<fastjet::contrib::SoftDrop>().delta_R();
	    zg = sdjets[i].structure_of<fastjet::contrib::SoftDrop>().symmetry();
	    mg = sdjets[i].m();
	}	    
	    
        mytxtfile << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << xsecweight << ", " << missing << ", " << jets[i].perp() << ", " << jets[i].m() << ", " << jets[i].constituents().size() << ", " << jets[i].eta() << ", " << jets[i].phi() << ", " << parton << ", " << mg << ", " << rg << ", " << zg << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << -999 << ", " << i << ", " << jets[i].delta_R(partons.at(0)) <<", " << jets[i].delta_R(partons.at(1)) << ", " << theEvent->event_scale() << "\n";}
      
      for (unsigned int is = 0; is < jets_s.size(); is++)
      {
        if (jets_s[is].user_info<fastjet::MyUserInfo>().matching() == matched)
          continue;
	// smeared jets without true jets matched
        vector <double> jet_q_s = jets_s[is].user_info<fastjet::MyUserInfo>().q();
	int parton_s = jets_s[is].user_info<fastjet::MyUserInfo>().parton();
        jets_s[is].set_user_info(new fastjet::MyUserInfo(jet_q_s, fake, parton_s));
 	mytxtfile << "-999, -999, -999, -999, -999, -999, -999, -999, ";
 	for (vector<double>::iterator it = jet_q_s.begin(); it != jet_q_s.end(); ++it){
		    mytxtfile << *it << ", "; 
	}
	double rgs, zgs, mgs;
       	if (sdjets_s[is] == 0) // did not pass the softdrop requirement
	{
		    rgs = -2;
		    zgs = -2;
		    mgs = -2;
	}
	else {
		    rgs = sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().delta_R();
		    zgs = sdjets_s[is].structure_of<fastjet::contrib::SoftDrop>().symmetry();
		    mgs = sdjets_s[is].m();
	}	    
	    
	mytxtfile << xsecweight << ", " << fake << ", -999, -999, -999, -999, -999, -999, -999, -999, " << jets_s[is].perp() << ", " << jets_s[is].m() << ", " << jets_s[is].constituents().size() << ", " << jets_s[is].eta() << ", " << jets_s[is].phi() << ", " << parton_s << ", " << mgs << ", " << rgs << ", " << zgs << ", -999, -999, -999, -999\n";
      }*/

   }

    /// Normalise histograms etc., after the run
    void finalize()
    {
     mytxtfile.close();
    }

    ///@}

  private:
    double eff_array[45][80];
    std::ofstream mytxtfile;
    ifstream efftxtfile;
  };

  RIVET_DECLARE_PLUGIN(PYTHIA_JETS_CH);
}
