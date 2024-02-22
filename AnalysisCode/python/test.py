import pickle

def getMuonScaleFactors():
    fID = open('./../py_pkl/MuonEfficiencies_Run2012ReReco_53X.pkl','r')
    fISO = open('./../py_pkl/MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl','r')
    fTRIG = open('./../py_pkl/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.pkl','r')

    fDimuonsCorr = open("./../configs/DimuonCorr.txt", 'w')

    dataID = pickle.load(fID)['Loose']
    dataISO = pickle.load(fISO)['combRelIsoPF04dBeta<02_Loose']
    dataTRIG = pickle.load(fTRIG)['Mu17Mu8_OR_Mu17TkMu8']['Loose']['(eta,eta)']['(20<mu1<Infty,20<mu2<Infty)']

    ptBins = dataID['ptabseta1.2-2.1'].keys()
#    etaBins = ['ptabseta<0.9','ptabseta0.9-1.2','ptabseta1.2-2.1','ptabseta2.1-2.4']


    dataIDKey = dataID.keys()
#    print ('dataID key??')
#    print (len(dataIDKey))

    etaBins = dataID['abseta_2p4pt20-500'].keys()

    print ( 'MuonID - 0.9_1.2 : ',dataID['abseta_2p4pt20-500']['0.9_1.2']['data/mc']['efficiency_ratio'])
    print ( 'MuonID - 0.0,0.9 : ',dataID['abseta_2p4pt20-500']['0_0.9']['data/mc']['efficiency_ratio'])
    print ( 'MuonISO - 0.9_1.2 : ',dataISO['abseta_2p4pt20-500']['0.9_1.2']['data/mc']['efficiency_ratio'])
    print ( 'MuonISO - 0.0,0.9 : ',dataISO['abseta_2p4pt20-500']['0_0.9']['data/mc']['efficiency_ratio'])
    print ('DoubleMuonTRig 0.9,1.2 : ',dataTRIG['(0.9,1.2)(0.0,0.9)']['ratio']['efficiency'])
#    print (etaBins)
#    print ('finish')
#    print (dataTRIG)

#   print (len(dataTRIG))
    etamatching = { 
                   '0.9_1.2': '(0.9,1.2)',
                   '0_0.9'  : '(0.0,0.9)',
                   '1.2_2.1': '(1.2,2.1)', 
                   '2.1_2.4': '(2.1,2.4)'}

    dataTRIGKey = dataTRIG.keys()

    lep1idiso = [1 for _ in range(10)] 
    lep2idiso = [1 for _ in range(10)] 
    trigeff = [1 for _ in range(10)]
#    print (dataTRIGKey)
    for trigeta in range( len(dataTRIGKey) ) :
#        print(dataTRIGKey[trigeta])
#        print('00oo00oo')
        triggercorr = dataTRIG[dataTRIGKey[trigeta]]['ratio']['efficiency']
        trigeff[trigeta] = dataTRIG[dataTRIGKey[trigeta]]['ratio']['efficiency']
#        print (triggercorr)
           
        for etabin in range( len(etaBins) ) :
#            print ('etabin start')
#            print ( etaBins[etabin] )
            idcorr = dataID['abseta_2p4pt20-500'][ etaBins[etabin] ]['data/mc']['efficiency_ratio']
            isocorr = dataISO['abseta_2p4pt20-500'][ etaBins[etabin] ]['data/mc']['efficiency_ratio']
#            print ('abseta_2p4pt20-500 ? id ? : ',idcorr,'iso ? : ', isocorr )
#            print (etamatching[etaBins[etabin] ])
#            print ( 'ddd ',dataTRIGKey[trigeta].count( etamatching[etaBins[etabin] ] ) )
            lepidisocorr = 0
            if (dataTRIGKey[trigeta].find( etamatching[etaBins[etabin] ] ) == 0) :
#                print ( dataTRIGKey[trigeta].count( etamatching[etaBins[etabin] ] ) )
                print( 'idcorr ? : ',idcorr)
                lepidisocorr = pow(idcorr*isocorr,1)
#                lep1idiso[trigeta] = idcorr*isocorr
                lep1idiso[trigeta] = idcorr*1
            if (dataTRIGKey[trigeta].find( etamatching[etaBins[etabin] ],5 ) > 4) :
                print ('--???--??--??')
                lepidisocorr = idcorr*isocorr
#                lep2idiso[trigeta] = idcorr*isocorr
                lep2idiso[trigeta] = idcorr*1
            pass
        lepcorr = triggercorr*lepidisocorr
#        fDimuonsCorr.write(etarangetxt)
#        fDimuonsCorr.write(lepcorrtxt)
    
    for j in range(10) :
        print'DoubleMuonEta_Range_%s : %s ' % (j, dataTRIGKey[j] )
        etarangetxt = 'DoubleMuonEta_Range_%s : %s \n' % (trigeta, dataTRIGKey[j] ) 


        print ('lep 1 idiso ??  ', lep1idiso[j])
        print ('lep 2 idiso ??  ', lep2idiso[j])
        print ('trgieff ??      ', trigeff[j]  )


        print'DoubleMuonCorr_%s : %s ' % (j, trigeff[j]*lep1idiso[j]*lep2idiso[j])
        lepcorrtxt = 'DoubleMuonCorr_%s : %s \n' % (j, trigeff[j]*lep1idiso[j]*lep2idiso[j])
        fDimuonsCorr.write(etarangetxt)
        fDimuonsCorr.write(lepcorrtxt)

    fDimuonsCorr.close()
#    print (dataID.keys())
#print ('start')
getMuonScaleFactors()
