import gzip,csv,json
rows=list(csv.DictReader(gzip.open('/tmp/phen.tsv.bgz','rt',errors='replace'),delimiter='\t'))
idx={x['phenotype']:x for x in rows}
def nc(p):
    x=idx.get(p)
    if not x: return None
    try: return int(float(x['n_cases']))
    except: return None

FAM={
 'Asthma':[('20002_1111','self-report'),('22127','doctor-diagnosed Q'),('6152_8','touchscreen composite')],
 'Hayfever/rhinitis':[('20002_1387','self-report'),('6152_9','touchscreen composite')],
 'Hypertension':[('20002_1065','self-report'),('6150_4','touchscreen composite')],
 'Angina':[('20002_1074','self-report'),('I20','hospital ICD'),('I9_UAP','curated endpoint')],
 'Myocardial infarction':[('20002_1075','self-report'),('I21','hospital ICD'),('I9_MI','curated endpoint')],
 'Ischaemic heart disease':[('I25','hospital ICD'),('I9_IHD','curated endpoint'),('I9_CHD','curated endpoint'),('I9_CORATHER','curated endpoint')],
 'Stroke':[('20002_1081','self-report'),('I63','hospital ICD'),('I9_STR','curated endpoint'),('C_STROKE','curated endpoint')],
 'Atrial fibrillation':[('20002_1471','self-report'),('I48','hospital ICD'),('CARDIAC_ARRHYTM','curated endpoint')],
 'Venous thromboembolism':[('20002_1094','self-report'),('I80','hospital ICD'),('I9_VTE','curated endpoint'),('I9_DVTANDPULM','curated endpoint')],
 'Osteoarthritis':[('20002_1465','self-report'),('M17','hospital ICD'),('M16','hospital ICD'),('M13_ARTHROSIS','curated endpoint'),('KNEE_ARTHROSIS','curated endpoint')],
 'Back pain':[('20002_1294','self-report'),('M54','hospital ICD'),('M13_DORSALGIA','curated endpoint')],
 'Cataract':[('20002_1278','self-report'),('H26','hospital ICD'),('H25','hospital ICD'),('H7_LENS','curated endpoint')],
 'Cholelithiasis':[('20002_1162','self-report'),('K80','hospital ICD'),('K11_GALLBILPANC','curated endpoint')],
 'Gastro-oesophageal reflux':[('20002_1138','self-report'),('K21','hospital ICD'),('K20','hospital ICD')],
 'Diverticular disease':[('20002_1458','self-report'),('K57','hospital ICD')],
 'Appendicitis':[('20002_1502','self-report'),('K35','hospital ICD'),('K11_APPENDIX','curated endpoint')],
 'Kidney/ureter stone':[('20002_1197','self-report'),('N20','hospital ICD')],
 'Pneumonia':[('20002_1398','self-report'),('J18','hospital ICD'),('PNEUMONIA','curated endpoint')],
 'Rheumatoid arthritis':[('20002_1464','self-report'),('M13_RHEUMA','curated endpoint')],
 'Prostate cancer':[('C61','hospital ICD'),('C3_PROSTATE','curated endpoint')],
 'Breast cancer':[('C50','hospital ICD'),('C3_BREAST_3','curated endpoint')],
 'Non-melanoma skin cancer':[('C44','hospital ICD'),('C3_OTHER_SKIN','curated endpoint')],
 'Colon cancer':[('C18','hospital ICD'),('C3_COLON','curated endpoint')],
 'Carpal tunnel syndrome':[('G56','hospital ICD'),('G6_CARPTU','curated endpoint')],
 'Shoulder lesion':[('M75','hospital ICD'),('M13_SHOULDER','curated endpoint')],
 'Iron deficiency anaemia':[('D50','hospital ICD'),('D3_ANAEMIA_IRONDEF','curated endpoint')],
 'Sleep disorder':[('G47','hospital ICD'),('SLEEP','curated endpoint')],
 'Inguinal/other hernia':[('K40','hospital ICD'),('K11_HERNIA','curated endpoint')],
 'Retinal detachment':[('H33','hospital ICD'),('H7_RETINALDETACH','curated endpoint')],
 'Eyelid disorder':[('H02','hospital ICD'),('H7_EYELIDDIS','curated endpoint')],
 'Depression/mental disorder':[('20002_1286','self-report'),('KRA_PSY_ANYMENTAL','curated endpoint')],
 'Migraine/headache':[('20002_1265','self-report'),('R51','hospital ICD')],
 'Ulcerative colitis':[('K51','hospital ICD'),('COLITNONINFNAS','curated endpoint')],
 'Anaemia (other)':[('D64','hospital ICD'),('D3_OTHERANAEMIA','curated endpoint')],
 'Bronchitis/COPD':[('20002_1113','self-report'),('20002_1412','self-report'),('BRONCHITIS','curated endpoint')],
}
CONT={
 '135':'N self-reported illnesses (reporting propensity)',
 '137':'N medications taken',
 '2178':'Overall health rating',
 '2090':'Seen GP for nerves/anxiety/depression',
 'ICDMAIN_ANY_ENTRY':'Any hospital ICD episode (hospitalisation propensity)',
 'XXI_HEALTHFACTORS':'Health-service contact factors',
 '20127_irnt':'Neuroticism',
 '21001_irnt':'BMI',
 '189_irnt':'Townsend deprivation',
 '6138_1':'College/University degree',
 '20116_2':'Current smoker',
}
MINC=2000
out={}; drop=[]
for f,arms in FAM.items():
    keep=[]
    for p,lab in arms:
        c=nc(p)
        if p not in idx: drop.append((f,p,'absent')); continue
        if c is not None and c<MINC: drop.append((f,p,f'cases={c}')); continue
        keep.append({'pheno':p,'deftype':lab,'cases':c,'desc':idx[p]['description']})
    if len(keep)>=2: out[f]=keep
    else: drop.append((f,'FAMILY','<2 arms'))
allp=sorted({a['pheno'] for v in out.values() for a in v} | set(CONT))
json.dump({'families':out,'contaminants':CONT,'all':allp},open('/tmp/config.json','w'),indent=1)
print("families:",len(out),"arms:",sum(len(v) for v in out.values()),"contaminants:",len(CONT),"unique files:",len(allp))
print("dropped:",drop)
