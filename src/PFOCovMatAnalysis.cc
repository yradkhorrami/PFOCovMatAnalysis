#include "PFOCovMatAnalysis.h"
using namespace lcio ;
using namespace marlin ;
using namespace std ;

PFOCovMatAnalysis aPFOCovMatAnalysis;

PFOCovMatAnalysis::PFOCovMatAnalysis() :

Processor("PFOCovMatAnalysis"),
m_nRun(0),
m_nEvt(0),
m_nRunSum(0),
m_nEvtSum(0),
m_pTFile(NULL),
m_pTTree(NULL)
{

	_description = "This code makes plots for investigating impact of track refitting with true mass hypothesis + taking neautral Hadrons as massive PFOs in CovMat calculations";

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"inputPfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"ClusterMCTruthLinkCollection",
					"Name of input m_ClusterMCTruthLink Collection",
					m_ClusterMCTruthLinkCollection,
					std::string("ClusterMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthClusterLinkCollection",
					"Name of input MCTruthClusterLink Collection",
					m_MCTruthClusterLinkCollection,
					std::string("MCTruthClusterLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"TrackMCTruthLinkCollection",
					"Name of input TrackMCTruthLink Collection",
					m_TrackMCTruthLinkCollection,
					std::string("MarlinTrkTracksMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthTrackLinkCollection",
					"Name of input MCTruthTrackLink Collection",
					m_MCTruthTrackLinkCollection,
					std::string("MCTruthMarlinTrkTracksLink")
				);


	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("PFOCovMatAnalysis.root")
				);

	registerProcessorParameter(	"MinWeightTrackMCTruthLink" ,
					"Minimum acceptable weight for Track -> MCParticle Link"  ,
					m_MinWeightTrackMCTruthLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"MinWeightMCTruthTrackLink" ,
					"Minimum acceptable weight for MCParticle -> Track Link"  ,
					m_MinWeightMCTruthTrackLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"MinWeightClusterMCTruthLink" ,
					"Minimum acceptable weight for Cluster -> MCParticle Link"  ,
					m_MinWeightClusterMCTruthLink ,
					float(0.9f)
				);

	registerProcessorParameter(	"MinWeightMCTruthClusterLink" ,
					"Minimum acceptable weight for MCParticle -> Cluster Link"  ,
					m_MinWeightMCTruthClusterLink ,
					float(0.9f)
				);

}

void PFOCovMatAnalysis::init()
{
	streamlog_out(MESSAGE) << "   init called  " << std::endl;
	printParameters();
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	m_pTFile = new TFile( m_rootFile.c_str() , "recreate" );
	m_pTTree = new TTree( "PFOCovMat" , "PFOCovMat" );
	m_pTTree->SetDirectory( m_pTFile );
	m_pTTree->Branch( "Event" , &m_nEvt , "Event/I" );
	m_pTTree->Branch( "Run" , &m_nRun , "Run/I" );
	m_pTTree->Branch( "PDG" , &m_PDG );
	m_pTTree->Branch( "Type" , &m_Type );
	m_pTTree->Branch( "nTracks" , &m_nTracks );
	m_pTTree->Branch( "charge" , &m_charge );
	m_pTTree->Branch( "mcpPx" , &m_mcpPx );
	m_pTTree->Branch( "mcpPy" , &m_mcpPy );
	m_pTTree->Branch( "mcpPz" , &m_mcpPz );
	m_pTTree->Branch( "mcpPt" , &m_mcpPt );
	m_pTTree->Branch( "mcpEnergy" , &m_mcpEnergy );
	m_pTTree->Branch( "mcpTheta" , &m_mcpTheta );
	m_pTTree->Branch( "mcpPhi" , &m_mcpPhi );
	m_pTTree->Branch( "pfoPx" , &m_pfoPx );
	m_pTTree->Branch( "pfoPy" , &m_pfoPy );
	m_pTTree->Branch( "pfoPz" , &m_pfoPz );
	m_pTTree->Branch( "pfoPt" , &m_pfoPt );
	m_pTTree->Branch( "pfoEnergy" , &m_pfoEnergy );
	m_pTTree->Branch( "pfoTheta" , &m_pfoTheta );
	m_pTTree->Branch( "pfoPhi" , &m_pfoPhi );
	m_pTTree->Branch( "ResidualPx" , &m_ResidualPx );
	m_pTTree->Branch( "ResidualPy" , &m_ResidualPy );
	m_pTTree->Branch( "ResidualPz" , &m_ResidualPz );
	m_pTTree->Branch( "ResidualEnergy" , &m_ResidualEnergy );
	m_pTTree->Branch( "ResidualTheta" , &m_ResidualTheta );
	m_pTTree->Branch( "ResidualPhi" , &m_ResidualPhi );
	m_pTTree->Branch( "ErrorEnergy" , &m_ErrorEnergy );
	m_pTTree->Branch( "ErrorTheta" , &m_ErrorTheta );
	m_pTTree->Branch( "ErrorPhi" , &m_ErrorPhi );
	m_pTTree->Branch( "NormalizedResidualEnergy" , &m_NormalizedResidualEnergy );
	m_pTTree->Branch( "NormalizedResidualTheta" , &m_NormalizedResidualTheta );
	m_pTTree->Branch( "NormalizedResidualPhi" , &m_NormalizedResidualPhi );
	m_pTTree->Branch( "sigmaPx2" , &m_sigmaPx2 );
	m_pTTree->Branch( "sigmaPxPy" , &m_sigmaPxPy );
	m_pTTree->Branch( "sigmaPy2" , &m_sigmaPy2 );
	m_pTTree->Branch( "sigmaPxPz" , &m_sigmaPxPz );
	m_pTTree->Branch( "sigmaPyPz" , &m_sigmaPyPz );
	m_pTTree->Branch( "sigmaPz2" , &m_sigmaPz2 );
	m_pTTree->Branch( "sigmaPxE" , &m_sigmaPxE );
	m_pTTree->Branch( "sigmaPyE" , &m_sigmaPyE );
	m_pTTree->Branch( "sigmaPzE" , &m_sigmaPzE );
	m_pTTree->Branch( "sigmaE2" , &m_sigmaE2 );
	m_pTTree->Branch( "linkedMCPToTrack1_PDG" , &m_linkedMCPToTrack1_PDG );
	m_pTTree->Branch( "linkedMCPToTrack1_TrueMomentum" , &m_linkedMCPToTrack1_TrueMomentum );
	m_pTTree->Branch( "linkedMCPToTrack1_TrueEnergy" , &m_linkedMCPToTrack1_TrueEnergy );
	m_pTTree->Branch( "linkedMCPToTrack2_PDG" , &m_linkedMCPToTrack2_PDG );
	m_pTTree->Branch( "linkedMCPToTrack2_TrueMomentum" , &m_linkedMCPToTrack2_TrueMomentum );
	m_pTTree->Branch( "linkedMCPToTrack2_TrueEnergy" , &m_linkedMCPToTrack2_TrueEnergy );
}

void PFOCovMatAnalysis::Clear()
{
	m_PDG.clear();
	m_Type.clear();
	m_nTracks.clear();
	m_charge.clear();
	m_mcpPx.clear();
	m_mcpPy.clear();
	m_mcpPz.clear();
	m_mcpPt.clear();
	m_mcpEnergy.clear();
	m_mcpTheta.clear();
	m_mcpPhi.clear();
	m_pfoPx.clear();
	m_pfoPy.clear();
	m_pfoPz.clear();
	m_pfoPt.clear();
	m_pfoEnergy.clear();
	m_pfoTheta.clear();
	m_pfoPhi.clear();
	m_ResidualPx.clear();
	m_ResidualPy.clear();
	m_ResidualPz.clear();
	m_ResidualEnergy.clear();
	m_ResidualTheta.clear();
	m_ResidualPhi.clear();
	m_ErrorEnergy.clear();
	m_ErrorTheta.clear();
	m_ErrorPhi.clear();
	m_NormalizedResidualEnergy.clear();
	m_NormalizedResidualTheta.clear();
	m_NormalizedResidualPhi.clear();
	m_sigmaPx2.clear();
	m_sigmaPxPy.clear();
	m_sigmaPy2.clear();
	m_sigmaPxPz.clear();
	m_sigmaPyPz.clear();
	m_sigmaPz2.clear();
	m_sigmaPxE.clear();
	m_sigmaPyE.clear();
	m_sigmaPzE.clear();
	m_sigmaE2.clear();
	m_linkedMCPToTrack1_PDG.clear();
	m_linkedMCPToTrack1_TrueMomentum.clear();
	m_linkedMCPToTrack1_TrueEnergy.clear();
	m_linkedMCPToTrack2_PDG.clear();
	m_linkedMCPToTrack2_TrueMomentum.clear();
	m_linkedMCPToTrack2_TrueEnergy.clear();
}

void PFOCovMatAnalysis::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void PFOCovMatAnalysis::processEvent( EVENT::LCEvent *pLCEvent )
{
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;

	LCCollection *inputPfoCollection{};
	this->Clear();
	try
	{
		streamlog_out(MESSAGE) << "	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
		streamlog_out(MESSAGE) << "	///////////////////////////////////////////////// processing event " << m_nEvt << " in run " << m_nRun << " /////////////////////////////////////////////////////////////////////" << std::endl ;
		streamlog_out(MESSAGE) << "	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
		inputPfoCollection = pLCEvent->getCollection( m_inputPfoCollection );
		for (int i_pfo = 0; i_pfo < inputPfoCollection->getNumberOfElements() ; ++i_pfo)
		{
			ReconstructedParticle* pfo = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_pfo ) );
			MCParticle *linkedMCP = NULL;
			MCParticle* linkedMCPsToTrack1 = NULL;
			MCParticle* linkedMCPsToTrack2 = NULL;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	Checking PFO[ " << i_pfo << " ] with " << pfo->getTracks().size() << " Tracks and " << pfo->getClusters().size() << " Clusters , Charged = " << pfo->getCharge() << std::endl;
			if ( pfo->getTracks().size() == 0 )
			{
				linkedMCP = getMCParticleLinkedToCluster( pLCEvent , pfo->getClusters()[ 0 ] );
			}
			else if ( pfo->getTracks().size() == 1 )
			{
				linkedMCP = getMCParticleLinkedToTrack( pLCEvent , pfo->getTracks()[ 0 ] );
			}
			else if ( pfo->getTracks().size() == 2 )
			{
				if ( fabs( pfo->getCharge() ) < 0.5 )
				{
					MCParticle* linkedMCPsToTrack = NULL;
					//MCParticle* linkedMCPsToCluster = NULL;
					MCParticle* track1ParentMCP = NULL;
					MCParticle* track2ParentMCP = NULL;
					for ( unsigned int i_trk = 0 ; i_trk < pfo->getTracks().size() ; ++i_trk )
					{
						linkedMCPsToTrack = NULL;
						streamlog_out(DEBUG4) << "		Checking Track [ " << i_trk << " ] " << std::endl;
						linkedMCPsToTrack = getMCParticleLinkedToTrack( pLCEvent , pfo->getTracks()[ i_trk ] );
						if ( linkedMCPsToTrack != NULL )
						{
							if ( i_trk == 0 && linkedMCPsToTrack->getParents().size() != 0 ) track1ParentMCP = linkedMCPsToTrack->getParents()[ 0 ];
							if ( i_trk == 1 && linkedMCPsToTrack->getParents().size() != 0 ) track2ParentMCP = linkedMCPsToTrack->getParents()[ 0 ];
							if ( i_trk == 0 ) linkedMCPsToTrack1 = linkedMCPsToTrack;
							if ( i_trk == 1 ) linkedMCPsToTrack2 = linkedMCPsToTrack;
							streamlog_out(DEBUG4) << "		Track [ " << i_trk << " ] linked to a MCParticle: PDG = " << linkedMCPsToTrack->getPDG() << std::endl;
						}
						//for ( unsigned int i_clu = 0 ; i_clu < pfo->getClusters().size() ; ++i_clu )
						//{
						//	linkedMCPsToCluster = getMCParticleLinkedToCluster( pLCEvent , pfo->getClusters()[ i_clu ] );
						//	if ( linkedMCPsToTrack == linkedMCPsToCluster ) linkedMCP = linkedMCPsToTrack->getParents()[ 0 ];
						//}
					}
					if ( track1ParentMCP != NULL && track2ParentMCP != NULL )
					{
						streamlog_out(DEBUG4) << "		Both parents of linked MCPs are found" << std::endl;
						if ( track1ParentMCP == track2ParentMCP )
						{
							linkedMCP = track1ParentMCP;
						}
					}
					else
					{
						streamlog_out(DEBUG4) << "		At-least one of parents of linked MCPs are not found" << std::endl;
					}
				}
				else
				{
					Track* track1 = pfo->getTracks()[ 0 ];
					TrackState *track1StateAtFirstHit = track1->getTrackStates()[ 1 ];
					Track* track2 = pfo->getTracks()[ 1 ];
					TrackState *track2StateAtFirstHit = track2->getTrackStates()[ 1 ];
					if ( sqrt( pow( track1StateAtFirstHit->getReferencePoint()[ 0 ] , 2 ) + pow( track1StateAtFirstHit->getReferencePoint()[ 1 ] , 2 ) + pow( track1StateAtFirstHit->getReferencePoint()[ 2 ] , 2 ) ) < sqrt( pow( track2StateAtFirstHit->getReferencePoint()[ 0 ] , 2 ) + pow( track2StateAtFirstHit->getReferencePoint()[ 1 ] , 2 ) + pow( track2StateAtFirstHit->getReferencePoint()[ 2 ] , 2 ) ) )
					{
						linkedMCP = getMCParticleLinkedToTrack( pLCEvent , track1 );
					}
					else
					{
						linkedMCP = getMCParticleLinkedToTrack( pLCEvent , track2 );
					}
				}
			}
			if ( linkedMCP != NULL ) getNormalizedResiduals( pfo , linkedMCP , linkedMCPsToTrack1 , linkedMCPsToTrack2 );
		}
		m_pTTree->Fill();
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
	}
}

EVENT::MCParticle* PFOCovMatAnalysis::getMCParticleLinkedToTrack( EVENT::LCEvent *pLCEvent , EVENT::Track* inputTrack )
{
	LCRelationNavigator navTrackMCTruth( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
	LCRelationNavigator navMCTruthTrack( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );

	const EVENT::LCObjectVec& mcpvec = navTrackMCTruth.getRelatedToObjects( inputTrack );
	const EVENT::FloatVec&  mcpweightvec = navTrackMCTruth.getRelatedToWeights( inputTrack );
	MCParticle *linkedMCP = NULL;
	double maxweightTrackToMCP = 0.0;
	int iTrackToMCPmax = -1;
	int iMCPToTrackmax = -1;
	streamlog_out(DEBUG0) << "	Checking " << mcpvec.size() << " MCPs linked to Track" << std::endl;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		streamlog_out(DEBUG0) << "	Checking MCP[ " << i_mcp << " ] , MCP PDG = " << testMCP->getPDG() << " , link weight = " << mcp_weight << std::endl;
		if ( mcp_weight > maxweightTrackToMCP && mcp_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightTrackToMCP = mcp_weight;
			iTrackToMCPmax = i_mcp;
		}
	}
	if ( iTrackToMCPmax != -1 )
	{
		MCParticle *testMCP = (MCParticle *) mcpvec.at( iTrackToMCPmax );
		Track *linkedTrack;
		const EVENT::LCObjectVec& trackvec = navMCTruthTrack.getRelatedToObjects( testMCP );
		const EVENT::FloatVec&  trackweightvec = navMCTruthTrack.getRelatedToWeights( testMCP );
		streamlog_out(DEBUG0) << "	Found linked MCP, MCP PDG: " << testMCP->getPDG() << " , link weight = " << maxweightTrackToMCP << " , looking for " << trackvec.size() << " Track(s) linked back to this MCParticle" << std::endl;
		double maxweightMCPToTrack = 0.0;
		for ( unsigned int i_track = 0; i_track < trackvec.size(); i_track++ )
		{
			double track_weight = trackweightvec.at( i_track );
			streamlog_out(DEBUG0) << "	Checking track[ " << i_track << " ] , link weight = " << track_weight << std::endl;
			if ( track_weight > maxweightMCPToTrack  && track_weight >= m_MinWeightMCTruthTrackLink )
			{
				maxweightMCPToTrack = track_weight;
				iMCPToTrackmax = i_track;
			}
		}
		if ( iMCPToTrackmax != -1 )
		{
			linkedTrack = (Track *) trackvec.at( iMCPToTrackmax );
			if ( linkedTrack == inputTrack )
			{
				linkedMCP = testMCP;
				streamlog_out(DEBUG1) << "	Found a MCParticle (PDGCode: " << linkedMCP->getPDG() << ") 	linked to track " << std::endl;
				streamlog_out(DEBUG1) << "	Track_MCParticle link weight = " << maxweightTrackToMCP << " 	; 	MCParticle_Track link weight = " << maxweightMCPToTrack << std::endl;
			}
		}
	}
	return linkedMCP;
}

EVENT::MCParticle* PFOCovMatAnalysis::getMCParticleLinkedToCluster( EVENT::LCEvent *pLCEvent , EVENT::Cluster* inputCluster )
{
	LCRelationNavigator navClusterMCTruth( pLCEvent->getCollection( m_ClusterMCTruthLinkCollection ) );
	LCRelationNavigator navMCTruthCluster( pLCEvent->getCollection( m_MCTruthClusterLinkCollection ) );

	const EVENT::LCObjectVec& mcpvec = navClusterMCTruth.getRelatedToObjects( inputCluster );
	const EVENT::FloatVec&  mcpweightvec = navClusterMCTruth.getRelatedToWeights( inputCluster );
	MCParticle *linkedMCP = NULL;
	double maxweightClusterToMCP = 0.0;
	int iClusterToMCPmax = -1;
	int iMCPToClustermax = -1;
	streamlog_out(DEBUG0) << "	Checking " << mcpvec.size() << " MCPs linked to Cluster" << std::endl;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double mcp_weight = mcpweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		streamlog_out(DEBUG0) << "	Checking MCP[ " << i_mcp << " ] , MCP PDG = " << testMCP->getPDG() << " , link weight = " << mcp_weight << std::endl;
		if ( mcp_weight > maxweightClusterToMCP && mcp_weight >= m_MinWeightClusterMCTruthLink && fabs( testMCP->getCharge() ) < 0.5 )
		{
			maxweightClusterToMCP = mcp_weight;
			iClusterToMCPmax = i_mcp;
		}
	}
	if ( iClusterToMCPmax != -1 )
	{
		linkedMCP = (MCParticle *) mcpvec.at( iClusterToMCPmax );
		Cluster *linkedCluster;
		const EVENT::LCObjectVec& clustervec = navMCTruthCluster.getRelatedToObjects( linkedMCP );
		const EVENT::FloatVec&  clusterweightvec = navMCTruthCluster.getRelatedToWeights( linkedMCP );
		streamlog_out(DEBUG0) << "	Found linked MCP, MCP PDG: " << linkedMCP->getPDG() << " , link weight = " << maxweightClusterToMCP << " , looking for " << clustervec.size() << " Cluster(s) linked back to this MCParticle" << std::endl;
		double maxweightMCPToCluster = 0.;
		for ( unsigned int i_cluster = 0; i_cluster < clustervec.size(); i_cluster++ )
		{
			double cluster_weight = clusterweightvec.at( i_cluster );
			streamlog_out(DEBUG0) << "	Checking cluster[ " << i_cluster << " ] , link weight = " << cluster_weight << std::endl;
			if ( cluster_weight > maxweightMCPToCluster  && cluster_weight >= m_MinWeightMCTruthClusterLink )
			{
				maxweightMCPToCluster = cluster_weight;
				iMCPToClustermax = i_cluster;
			}
		}
		if ( iMCPToClustermax != -1 )
		{
			linkedCluster = (Cluster *) clustervec.at( iMCPToClustermax );
			if ( linkedCluster == inputCluster )
			{
				streamlog_out(DEBUG1) << "	Found a MCParticle (PDGCode: " << linkedMCP->getPDG() << ") 	linked to cluster " << std::endl;
				streamlog_out(DEBUG1) << "	Cluster_MCParticle link weight = " << maxweightClusterToMCP << " 	; 	MCParticle_Cluster link weight = " << maxweightMCPToCluster << std::endl;
			}
		}
	}
	return linkedMCP;
}

void PFOCovMatAnalysis::getNormalizedResiduals( EVENT::ReconstructedParticle* pfo , EVENT::MCParticle* linkedMCP , EVENT::MCParticle* linkedMCPToTrack1 , EVENT::MCParticle* linkedMCPToTrack2 )
{
	TVector3 mcpMomentum( linkedMCP->getMomentum() );
	double mcpEnergy = linkedMCP->getEnergy();
	TVector3 pfoMomentum( pfo->getMomentum() );
	double pfoEnergy = pfo->getEnergy();
	std::vector<float> pfoCovMat = pfo->getCovMatrix();
	TVector3 mcpMomentumUnit = mcpMomentum; mcpMomentumUnit.SetMag( 1.0 );
	TVector3 mcpPtUnit( mcpMomentum.Px() , mcpMomentum.Py() , 0.0 ); mcpPtUnit.SetMag( 1.0 );
	TVector3 pfoMomentumUnitRotated = pfoMomentum; pfoMomentumUnitRotated.SetPhi( mcpMomentum.Phi() ); pfoMomentumUnitRotated.SetMag( 1.0 );
	TVector3 pfoPtUnit( pfoMomentum.Px() , pfoMomentum.Py() , 0.0 ); pfoPtUnit.SetMag( 1.0 );

	double PxResidual	= ( pfoMomentum - mcpMomentum ).Px();
	double PyResidual	= ( pfoMomentum - mcpMomentum ).Py();
	double PzResidual	= ( pfoMomentum - mcpMomentum ).Pz();
	double energyResidual	= pfoEnergy - mcpEnergy;
	//double thetaResidual	= ( pfoMomentum.Theta() >= mcpMomentum.Theta() ? acos( mcpMomentumUnit.Dot( pfoMomentumUnitRotated ) ) : -1.0 * acos( mcpMomentumUnit.Dot( pfoMomentumUnitRotated ) ) );
	double thetaResidual	= pfoMomentum.Theta() - mcpMomentum.Theta();
	//double phiResidual	= ( pfoMomentum.Phi() >= mcpMomentum.Phi() ? acos( mcpPtUnit.Dot( pfoPtUnit ) ) : -1.0 * acos( mcpPtUnit.Dot( pfoPtUnit ) ) );
	double phiResidual	= pfoMomentum.Phi() - mcpMomentum.Phi();
	if ( M_PI < phiResidual ) phiResidual = phiResidual - 2.0 * M_PI;
	if ( phiResidual < -1.0 * M_PI ) phiResidual = phiResidual + 2.0 * M_PI;

	double Px		= pfoMomentum.Px();
	double Py		= pfoMomentum.Py();
	double Pz		= pfoMomentum.Pz();
	double P2		= pfoMomentum.Mag2();
	double Pt2		= pow( Px , 2 ) + pow( Py , 2 );
	double Pt		= sqrt( Pt2 );
	double sigmaPx2		= pfoCovMat[ 0 ];
	double sigmaPxPy	= pfoCovMat[ 1 ];
	double sigmaPy2		= pfoCovMat[ 2 ];
	double sigmaPxPz	= pfoCovMat[ 3 ];
	double sigmaPyPz	= pfoCovMat[ 4 ];
	double sigmaPz2		= pfoCovMat[ 5 ];
	double sigmaPxE		= pfoCovMat[ 6 ];
	double sigmaPyE		= pfoCovMat[ 7 ];
	double sigmaPzE		= pfoCovMat[ 8 ];
	double sigmaE2		= pfoCovMat[ 9 ];

	double dTheta_dPx	= Px * Pz / ( P2 * Pt );
	double dTheta_dPy	= Py * Pz / ( P2 * Pt );
	double dTheta_dPz	= -Pt / P2;
	double dPhi_dPx		= -Py / Pt2;
	double dPhi_dPy		= Px / Pt2;

	double sigmaE		= std::sqrt( sigmaE2 );
	double sigmaTheta	= std::sqrt( std::fabs( std::pow( dTheta_dPx , 2 ) * sigmaPx2 + std::pow( dTheta_dPy , 2 ) * sigmaPy2 + std::pow( dTheta_dPz , 2 ) * sigmaPz2 +
	 			2.0 * dTheta_dPx * dTheta_dPy * sigmaPxPy + 2.0 * dTheta_dPx * dTheta_dPz * sigmaPxPz + 2.0 * dTheta_dPy * dTheta_dPz * sigmaPyPz ) );
	double sigmaPhi		= std::sqrt( std::fabs( std::pow( dPhi_dPx , 2 ) * sigmaPx2 + std::pow( dPhi_dPy , 2 ) * sigmaPy2 + 2.0 * dPhi_dPx * dPhi_dPy * sigmaPxPy ) );

	m_PDG.push_back( linkedMCP->getPDG() );
	m_Type.push_back( pfo->getType() );
	m_nTracks.push_back( pfo->getTracks().size() );
	m_charge.push_back( pfo->getCharge() );
	m_mcpPx.push_back( linkedMCP->getMomentum()[ 0 ] );
	m_mcpPy.push_back( linkedMCP->getMomentum()[ 1 ] );
	m_mcpPz.push_back( linkedMCP->getMomentum()[ 2 ] );
	m_mcpPt.push_back( sqrt( pow( linkedMCP->getMomentum()[ 0 ] , 2 ) + pow( linkedMCP->getMomentum()[ 1 ] , 2 ) ) );
	m_mcpEnergy.push_back( linkedMCP->getEnergy() );
	m_mcpTheta.push_back( mcpMomentum.Theta() );
	m_mcpPhi.push_back( mcpMomentum.Phi() );
	m_pfoPx.push_back( pfo->getMomentum()[ 0 ] );
	m_pfoPy.push_back( pfo->getMomentum()[ 1 ] );
	m_pfoPz.push_back( pfo->getMomentum()[ 2 ] );
	m_pfoPt.push_back( sqrt( pow( pfo->getMomentum()[ 0 ] , 2 ) + pow( pfo->getMomentum()[ 1 ] , 2 ) ) );
	m_pfoEnergy.push_back( pfo->getEnergy() );
	m_pfoTheta.push_back( pfoMomentum.Theta() );
	m_pfoPhi.push_back( pfoMomentum.Phi() );
	m_ResidualPx.push_back( PxResidual );
	m_ResidualPy.push_back( PyResidual );
	m_ResidualPz.push_back( PzResidual );
	m_ResidualEnergy.push_back( energyResidual );
	m_ResidualTheta.push_back( thetaResidual );
	m_ResidualPhi.push_back( phiResidual );
	m_ErrorEnergy.push_back( sigmaE );
	m_ErrorTheta.push_back( sigmaTheta );
	m_ErrorPhi.push_back( sigmaPhi );
	m_NormalizedResidualEnergy.push_back( energyResidual / sigmaE );
	m_NormalizedResidualTheta.push_back( thetaResidual / sigmaTheta );
	m_NormalizedResidualPhi.push_back( phiResidual / sigmaPhi );
	m_sigmaPx2.push_back( sigmaPx2 );
	m_sigmaPxPy.push_back( sigmaPxPy );
	m_sigmaPy2.push_back( sigmaPy2 );
	m_sigmaPxPz.push_back( sigmaPxPz );
	m_sigmaPyPz.push_back( sigmaPyPz );
	m_sigmaPz2.push_back( sigmaPz2 );
	m_sigmaPxE.push_back( sigmaPxE );
	m_sigmaPyE.push_back( sigmaPyE );
	m_sigmaPzE.push_back( sigmaPzE );
	m_sigmaE2.push_back( sigmaE2 );
	if ( linkedMCPToTrack1 != NULL && linkedMCPToTrack2 != NULL  )
	{
		m_linkedMCPToTrack1_PDG.push_back( linkedMCPToTrack1->getPDG() );
		m_linkedMCPToTrack1_TrueMomentum.push_back( sqrt( pow( linkedMCPToTrack1->getMomentum()[ 0 ] , 2 ) + pow( linkedMCPToTrack1->getMomentum()[ 0 ] , 1 ) + pow( linkedMCPToTrack1->getMomentum()[ 2 ] , 2 ) ) );
		m_linkedMCPToTrack1_TrueEnergy.push_back( linkedMCPToTrack1->getEnergy() );
		m_linkedMCPToTrack2_PDG.push_back( linkedMCPToTrack2->getPDG() );
		m_linkedMCPToTrack2_TrueMomentum.push_back( sqrt( pow( linkedMCPToTrack2->getMomentum()[ 0 ] , 2 ) + pow( linkedMCPToTrack2->getMomentum()[ 0 ] , 1 ) + pow( linkedMCPToTrack2->getMomentum()[ 2 ] , 2 ) ) );
		m_linkedMCPToTrack2_TrueEnergy.push_back( linkedMCPToTrack2->getEnergy() );
	}
	else
	{
		m_linkedMCPToTrack1_PDG.push_back( 0 );
		m_linkedMCPToTrack1_TrueMomentum.push_back( 0.0 );
		m_linkedMCPToTrack1_TrueEnergy.push_back( 0.0 );
		m_linkedMCPToTrack2_PDG.push_back( 0 );
		m_linkedMCPToTrack2_TrueMomentum.push_back( 0.0 );
		m_linkedMCPToTrack2_TrueEnergy.push_back( 0.0 );
	}
}

void PFOCovMatAnalysis::check()
{

}

void PFOCovMatAnalysis::end()
{
	m_pTFile->cd();
	m_pTTree->Write();
	m_pTFile->Close();
	delete m_pTFile;
}
