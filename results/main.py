import pandas as pd

def process_csv(input_file, instance_group, method):
    data = pd.read_csv(input_file, header=None)
    data = data.iloc[:, -5:]  # Select the last 5 columns
    data.columns = ['instance', 'alfa', 'obj', 'timlim', 'time']

    data_cv = data.groupby(['instance', 'alfa'])['obj'].agg(['mean', 'std', 'min']).reset_index()
    data_cv['cv'] = data_cv['std'] / data_cv['mean']
    data_cv['instance_prefix'] = data_cv['instance'].astype(str).str[:2]

    data_cv.groupby(['instance_prefix', 'alfa'], as_index=False)[['cv', 'min']].mean() #.to_csv(output_file, index=False)
    data_cv['method'] = method
    data_cv['instance_group'] = instance_group
    data_cv = data_cv.fillna(0)
    return data_cv

# Lista de arquivos e parâmetros
files = [
    ('Projetos/FJS/complete-ILS-CN-NewRLS.csv', 'complete', 'ILS-CN'),
    ('Projetos/FJS/complete-ILS-RN-NewRLS.csv', 'complete', 'ILS-RN'),
    ('Projetos/FJS/medium-TS-CN-NewLS.csv', 'complete', 'TS-CN'),
    ('Projetos/FJS/medium-TS-RN-NewLS.csv', 'complete', 'TS-RN'),
    ('Projetos/FJS/complete-grasp-CN-definitivo.csv', 'complete', 'grasp-CN'),
    ('Projetos/FJS/complete-grasp-RN-definitivo.csv', 'complete', 'grasp-RN'),
    ('Projetos/FJS/medium-SA-fix.csv', 'complete', 'SA'),
    ('Projetos/FJS/mini-ILS-CN-definitivo.csv', 'mini', 'ILS-CN'),
    ('Projetos/FJS/mini-ILS-RN-definitivo.csv', 'mini', 'ILS-RN'),
    ('Projetos/FJS/mini-TS-CN-definitivo.csv', 'mini', 'TS-CN'),
    ('Projetos/FJS/mini-TS-RN-definitivo.csv', 'mini', 'TS-RN'),
    ('Projetos/FJS/mini-grasp-CN-definitivo.csv', 'mini', 'grasp-CN'),
    ('Projetos/FJS/mini-grasp-RN-definitivo.csv', 'mini', 'grasp-RN'),
    ('Projetos/FJS/mini-SA-definitivo.csv', 'mini', 'SA'),
    ('Projetos/FJS/classic-ILS-CN-definitivo.csv', 'classic', 'ILS-CN'),
    ('Projetos/FJS/classic-ILS-RN-definitivo.csv', 'classic', 'ILS-RN'),
    ('Projetos/FJS/classic-TS-CN-definitivo.csv', 'classic', 'TS-CN'),
    ('Projetos/FJS/classic-TS-RN-definitivo.csv', 'classic', 'TS-RN'),
    ('Projetos/FJS/classic-grasp-CN-definitivo.csv', 'classic', 'grasp-CN'),
    ('Projetos/FJS/classic-grasp-RN-definitivo.csv', 'classic', 'grasp-RN'),
    ('Projetos/FJS/classic-SA-definitivo.csv', 'classic', 'SA'),
    ('Projetos/FJS/classic-Tayebi.csv', 'classic', 'Tayebi')
]

# Concatenando os DataFrames
df_combined = pd.concat([process_csv(f, c, m) for f, c, m in files], ignore_index=True)
df_combined.to_csv("results.csv", index=False)