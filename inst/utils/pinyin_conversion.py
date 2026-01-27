import pyreadr
from pypinyin import lazy_pinyin
#print(''.join(lazy_pinyin('地龙')))

# load data
data = pyreadr.read_r('./herb_data.rda')
result = data['herb_data']
result['Herb_cn_name'] = result['Herb_cn_name'].astype(str)

# use pypinyin to convert Chinese name into pinyin name
result['Herb_pinyin_name'] = result['Herb_cn_name'].apply(lambda x: ''.join(lazy_pinyin(x)))

# output
pyreadr.write_rdata('./herb_data.rda', result)