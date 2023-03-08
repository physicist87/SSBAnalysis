import pickle

def getMuonScaleFactors():
    fID = open('MuonEfficiencies_Run2012ReReco_53X.pkl','r')
    fISO = open('MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl','r')
    fTRIG = open('MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.pkl','r')

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
#    print (dataTRIGKey)
    for trigeta in range( len(dataTRIGKey) ) :
#        print(dataTRIGKey[trigeta])
#        print('00oo00oo')
        triggercorr = dataTRIG[dataTRIGKey[trigeta]]['ratio']['efficiency']
#        print (triggercorr)
           
        for etabin in range( len(etaBins) ) :
#            print ('etabin start')
#            print ( etaBins[etabin] )
            idcorr = dataID['abseta_2p4pt20-500'][ etaBins[etabin] ]['data/mc']['efficiency_ratio']
            isocorr = dataID['abseta_2p4pt20-500'][ etaBins[etabin] ]['data/mc']['efficiency_ratio']
#            print ('abseta_2p4pt20-500 ? id ? : ',idcorr,'iso ? : ', isocorr )
#            print (etamatching[etaBins[etabin] ])
#            print ( 'ddd ',dataTRIGKey[trigeta].count( etamatching[etaBins[etabin] ] ) )
            lepidisocorr = 1
            if (dataTRIGKey[trigeta].count( etamatching[etaBins[etabin] ] ) == 2) :
#                print ( dataTRIGKey[trigeta].count( etamatching[etaBins[etabin] ] ) )
                lepidisocorr = pow(idcorr*idcorr,2)
            elif (dataTRIGKey[trigeta].count( etamatching[etaBins[etabin] ] ) == 1) :
#                print ('--???--??--??')
                lepidisocorr = pow(idcorr*idcorr,1)
            else :
                lepidisocorr =1
                pass
            pass
        lepcorr = triggercorr*lepidisocorr
        print'DoubleMuonEta_Range_%s : %s ' % (trigeta, dataTRIGKey[trigeta] )
        etarangetxt = 'DoubleMuonEta_Range_%s : %s \n' % (trigeta, dataTRIGKey[trigeta] ) 
        print'DoubleMuonCorr_%s : %s ' % (trigeta, lepcorr)
        lepcorrtxt = 'DoubleMuonCorr_%s : %s \n' % (trigeta, lepcorr)
        fDimuonsCorr.write(etarangetxt)
        fDimuonsCorr.write(lepcorrtxt)
    fDimuonsCorr.close()
#    print (dataID.keys())
#print ('start')
getMuonScaleFactors()
