# A03.flowscale
1. How to get pt integrated vn from the masured dN/dpt and vn(pt)
2. How to get eta integrated vn from the measured dvn/d\eta 

inputs : 

-vn measurements  s_nn = 5.02TeV
https://www.hepdata.net/record/ins1666817
 v2{2,|Δη|>1.} as a function of pT

-pt spectra s_nn = 5.02TeV
https://www.hepdata.net/record/86210
#frac{1}{N_{evt}} #frac{d^{2} N}{dp_{T}d#eta}(Gev^{-1}c) 

-vn(\eta) s_nn = 2.76TeV, getting 5.02TeV data from Freja soon
https://www.hepdata.net/record/ins1456145

Exercise 1. 
Calculate  pt-integrated v2 with the masured dN/dpt and vn(pt) as a function of the centrality.
You can find the necessary inputs from the following HEPdata for \sqrt{s_{NN}} = 5.02TeV
(1)Input dN/dpt : pt spectra  https://www.hepdata.net/record/86210 
(2)Input v2{2,|Δη|>1.} as a function of pT https://www.hepdata.net/record/ins1666817

Q1. Calculate the pt-integrated v2 both for (a)0.2<p_{T}<3.0GeV/c and (b)0.<p_{T}<3.0GeV/c and compare to (c) Table 1 from (2).
Q2. Draw (a), (b) and (c) and draw the ratios (a)/(c) and (b)/(c).

Exercise 2. 
v2{2}(\eta) data can be found in https://doi.org/10.17182/hepdata.73940 for various centrality ranges. 
Please calculate the \eta-integrated v2{2} for three ALICE detectors as a function of centrality.
Here is the \eta acceptance of three detectors,  TPC: {-0.8,0.8}, V0A:{2.8,5.1},  V0C:{-3.7,-1.7} 
Q1. Draw the \eta-integrated v2{2} for three ALICE detectors as a function of centrality.
Q2. Draw the ratios to TPC: {-0.8,0.8}
