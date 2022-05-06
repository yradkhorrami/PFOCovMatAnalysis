#ifndef PFOCovMatAnalysis_h
#define PFOCovMatAnalysis_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/Cluster.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;
class TF1;

using namespace lcio ;
using namespace marlin ;

class PFOCovMatAnalysis : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new PFOCovMatAnalysis;
		}
		PFOCovMatAnalysis();
		virtual ~PFOCovMatAnalysis() = default;
		PFOCovMatAnalysis(const PFOCovMatAnalysis&) = delete;
		PFOCovMatAnalysis& operator=(const PFOCovMatAnalysis&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		EVENT::MCParticle* getMCParticleLinkedToTrack( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrack );
		EVENT::MCParticle* getMCParticleLinkedToCluster( EVENT::LCEvent *pLCEvent , EVENT::Cluster* inputCluster );
		virtual void getNormalizedResiduals( EVENT::ReconstructedParticle* pfo , EVENT::MCParticle* linkedMCP );
		virtual void check();
		virtual void end();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_inputPfoCollection{};
		std::string				m_TrackMCTruthLinkCollection{};
		std::string				m_MCTruthTrackLinkCollection{};
		std::string				m_ClusterMCTruthLinkCollection{};
		std::string				m_MCTruthClusterLinkCollection{};
		std::string				m_rootFile{};

		float					m_MinWeightTrackMCTruthLink{};
		float					m_MinWeightMCTruthTrackLink{};
		float					m_MinWeightClusterMCTruthLink{};
		float					m_MinWeightMCTruthClusterLink{};

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		TFile					*m_pTFile;
	        TTree					*m_pTTree;

		std::vector<int>			m_PDG{};
		std::vector<int>			m_Type{};
		std::vector<int>			m_nTracks{};
		std::vector<float>			m_charge{};
		std::vector<float>			m_mcpPx{};
		std::vector<float>			m_mcpPy{};
		std::vector<float>			m_mcpPz{};
		std::vector<float>			m_mcpPt{};
		std::vector<float>			m_mcpEnergy{};
		std::vector<float>			m_mcpTheta{};
		std::vector<float>			m_mcpPhi{};
		std::vector<float>			m_pfoPx{};
		std::vector<float>			m_pfoPy{};
		std::vector<float>			m_pfoPz{};
		std::vector<float>			m_pfoPt{};
		std::vector<float>			m_pfoEnergy{};
		std::vector<float>			m_pfoTheta{};
		std::vector<float>			m_pfoPhi{};
		std::vector<float>			m_ResidualEnergy{};
		std::vector<float>			m_ResidualTheta{};
		std::vector<float>			m_ResidualPhi{};
		std::vector<float>			m_ErrorEnergy{};
		std::vector<float>			m_ErrorTheta{};
		std::vector<float>			m_ErrorPhi{};
		std::vector<float>			m_NormalizedResidualEnergy{};
		std::vector<float>			m_NormalizedResidualTheta{};
		std::vector<float>			m_NormalizedResidualPhi{};


};
#endif
