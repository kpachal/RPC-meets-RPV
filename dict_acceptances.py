# These contain two dodgy extrapolations
# 1) Changing 8 TeV -> 13 TeV means different pdf effects,
#    so acceptances actually drop a lot faster at 8 TeV
#    at high mass. However, at low mass, difference is small
# 2) The mjj cut in 8 TeV dijets was lower than the
#    corresponding cut in TLA, so the acceptance is artifically
#    enlarged.
# I am assuming that these two effects roughly cancel out.

dict_acceptances = {
"500" : 0.48,
"600" : 0.46,
"700" : 0.45,
"800" : 0.43,
"900" : 0.42,
"1000" : 0.41,
"1100" : -1,
"1200" : 0.38,
"1300" : -1,
"1400" : -1,
"1500" : 0.41,
"1600" : -1,
"1700" : -1,
"1800" : 0.38,
"1900" : -1,
"2000" : 0.37,
"2100" : -1,
"2200" : -1,
"2300" : -1,
"2400" : 0.33,
"2500" : -1,
"2600" : -1,
"2700" : -1,
"2800" : 0.30,
"2900" : -1,
"3000" : -1,
"3100" : -1,
"3200" : 0.26,
"3600" : 0.21,
"4000" : 0.17
}
