import pickle
data1 = pickle.load(open("MuonEfficiencies_Run2012ReReco_53X.pkl","rb"))
data2 = pickle.load(open("MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl","rb"))
data3 = pickle.load(open("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl","rb"))
data4 = pickle.load(open("MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.pkl","rb"))
data5 = pickle.load(open("MuHLTEfficiencies_Run_2012ABCD_53X_DR01-2.pkl","rb"))

output1 = open("MuonEfficiencies_Run2012ReReco_53X.txt","w")
output2 = open("MuonEfficiencies_ISO_Run_2012ReReco_53X.txt","w")
output3 = open("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.txt","w")
output4 = open("MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.txt","w")
output5 = open("MuHLTEfficiencies_Run_2012ABCD_53X_DR01-2.txt","w")

output1.write(str(data1))
output2.write(str(data2))
output3.write(str(data3))
output4.write(str(data4))
output5.write(str(data5))

output1.flush()
output2.flush()
output3.flush()
output4.flush()
output5.flush()

output1.close()
output2.close()
output3.close()
output4.close()
output5.close()
