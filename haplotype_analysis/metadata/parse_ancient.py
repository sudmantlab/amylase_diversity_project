import pdb
import json
import pandas as pd
f = open("./config_final_ancient_info.json",'r')
d = json.load(f)


outrows = []
for datum in d.items():
    sample = datum[0]
    outrows.append({"sample":sample,
                    "popGrouping":datum[1]['popGrouping'],
                    "ageAverage":datum[1]['ageAverage']})
    
t = pd.DataFrame(outrows)
t.to_csv("ancient_info.tsv",sep="\t", index=False)
