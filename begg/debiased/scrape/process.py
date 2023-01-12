from lxml import etree
import glob
import pdb
import re
import copy
import json

out_file = 'cochrane_data.json'
save_path = './download/*.rm5'
filenames = glob.glob(save_path)

# j = 2
dataset = []

for filename in filenames:
    # file = './download/CD002258StatsDataOnly.rm5'
    # pdb.set_trace()
    # print(filename)
    root = etree.parse(filename).getroot()

    authors = [str(person.find('LAST_NAME').text) + ', ' + str(person.find('FIRST_NAME').text) 
               for person in list(root.find('COVER_SHEET').find('CREATORS')) if person.tag=='PERSON']
    # title = root.find('COVER_SHEET').find('TITLE').text
    # print(filename)
    # if filename=='./download/CD000031StatsDataOnly.rm5':
    #     pdb.set_trace()
    title = re.search('<TITLE[^>]*>(.*)</TITLE>',str(etree.tostring(root.find('COVER_SHEET').find('TITLE')))).group(1)
    # ID = re.search('(*[0-9]+)',file).group(1)
    ID = filename.split('/')[2][:-17]
    # if ID == "CD001423":
    #     pdb.set_trace()
    # outcomes = [{'n_studies' : int(outcome.attrib['STUDIES']),
    #              'effect_measure' : outcome.attrib['EFFECT_MEASURE'],
    #              'outcome_pval' : float(outcome.attrib['P_Z']),
    #              'outcome_name' : outcome.find('NAME').text,
    #              'outcome_id' : outcome.attrib['ID'],
    #              'comparison_name' : comparison.find('NAME').text,
    #              'comparison_id' : comparison.attrib['ID'],
    #              'studies' :  [{'CI' : [float(data.attrib['CI_START']),float(data.attrib['CI_END'])],
    #                             'effect' : float(data.attrib['EFFECT_SIZE']),
    #                             'estimable' : data.attrib['ESTIMABLE'],
    #                             'se' : float(data.attrib['SE']),
    #                             'name' : data.attrib['STUDY_ID']
    #              }
    #                            for data in outcome.xpath('.//DICH_DATA')+outcome.xpath('.//CONT_DATA')]
    # }
    #             for comparison in root.find('ANALYSES_AND_DATA')
    #             for outcome in comparison.findall('DICH_OUTCOME') + comparison.findall('CONT_OUTCOME')] # omitting IV outcomes
    # pdb.set_trace()
    outcomes = {}
    for comparison in root.find('ANALYSES_AND_DATA'):
        for outcome in comparison.findall('DICH_OUTCOME') + comparison.findall('CONT_OUTCOME'):
            studies = {}
            for data in outcome.xpath('.//DICH_DATA')+outcome.xpath('.//CONT_DATA'):
                study_ID = data.attrib['STUDY_ID']
                if not study_ID in studies:
                    studies[study_ID] =  {'CI' : [float(data.attrib['CI_START']),float(data.attrib['CI_END'])],
                                  'effect' : float(data.attrib['EFFECT_SIZE']),
                                  'estimable' : data.attrib['ESTIMABLE'],
                                  'se' : float(data.attrib['SE']),
                                  'ID' : data.attrib['STUDY_ID']}
            outcomes[outcome.attrib['ID']] = {'n_studies' : int(outcome.attrib['STUDIES']),
                          'effect_measure' : outcome.attrib['EFFECT_MEASURE'],
                          'outcome_pval' : float(outcome.attrib['P_Z']),
                          'outcome_name' : outcome.find('NAME').text,
                          'outcome_id' : outcome.attrib['ID'],
                          'comparison_name' : comparison.find('NAME').text,
                          'comparison_id' : comparison.attrib['ID'],
                          'studies' : studies}
    # pdb.set_trace()
    dataset.append({'authors' : authors, 'title' : title, 'ID' : ID, 'outcomes' : outcomes})



## only mean difference effects and n>1 metaanalyses
out = copy.deepcopy(dataset)
for paper in out:
    paper['outcomes'] = dict(filter(lambda item : item[1]['effect_measure']=='MD' and item[1]['n_studies'] > 1,paper['outcomes'].items()))
out = list(filter(lambda paper : len(paper['outcomes']) > 0, out))

with open(out_file,'w') as f:
    f.write(json.dumps(out))
