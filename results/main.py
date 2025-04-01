import pandas as pd

def process_csv(input_file, output_file):
    data = pd.read_csv(input_file, header=None)
    data = data.iloc[:, -5:]  # Select the last 5 columns
    data.columns = ['instance', 'alfa', 'obj', 'timlim', 'time']

    data_cv = data.groupby(['instance', 'alfa'])['obj'].agg(['mean', 'std', 'min']).reset_index()
    data_cv['cv'] = data_cv['std'] / data_cv['mean']
    data_cv['instance_prefix'] = data_cv['instance'].astype(str).str[:2]

    data_cv.groupby(['instance_prefix', 'alfa'], as_index=False)[['cv', 'min']].mean().to_csv(output_file, index=False)

# Example usage
process_csv('Projetos/FJS/complete-ILS-CN-NewRLS.csv', 'summary-complete-ILS-CN.csv')
process_csv('Projetos/FJS/complete-ILS-RN-NewRLS.csv', 'summary-complete-ILS-RN.csv')
process_csv('Projetos/FJS/complete-TS-CN-NewLS.csv', 'summary-complete-TS-CN.csv')
process_csv('Projetos/FJS/complete-TS-RN-NewLS.csv', 'summary-complete-TS-RN.csv')
process_csv('Projetos/FJS/complete-grasp-CN-definitivo.csv', 'summary-complete-grasp-CN.csv')
process_csv('Projetos/FJS/complete-grasp-RN-definitivo.csv', 'summary-complete-grasp-RN.csv')
process_csv('Projetos/FJS/medium-SA-fix.csv', 'summary-complete-SA.csv')

process_csv('Projetos/FJS/mini-ILS-CN-definitivo.csv', 'summary-mini-ILS-CN.csv')
process_csv('Projetos/FJS/mini-ILS-RN-definitivo.csv', 'summary-mini-ILS-RN.csv')
process_csv('Projetos/FJS/mini-TS-CN-definitivo.csv', 'summary-mini-TS-CN.csv')
process_csv('Projetos/FJS/mini-TS-RN-definitivo.csv', 'summary-mini-TS-RN.csv')
process_csv('Projetos/FJS/mini-grasp-CN-definitivo.csv', 'summary-mini-grasp-CN.csv')
process_csv('Projetos/FJS/mini-grasp-RN-definitivo.csv', 'summary-mini-grasp-RN.csv')
process_csv('Projetos/FJS/mini-SA-definitivo.csv', 'summary-mini-SA.csv')

process_csv('Projetos/FJS/classic-ILS-CN-definitivo.csv', 'summary-classic-ILS-CN.csv')
process_csv('Projetos/FJS/classic-ILS-RN-definitivo.csv', 'summary-classic-ILS-RN.csv')
process_csv('Projetos/FJS/classic-TS-CN-definitivo.csv', 'summary-classic-TS-CN.csv')
process_csv('Projetos/FJS/classic-TS-RN-definitivo.csv', 'summary-classic-TS-RN.csv')
process_csv('Projetos/FJS/classic-grasp-CN-definitivo.csv', 'summary-classic-grasp-CN.csv')
process_csv('Projetos/FJS/classic-grasp-RN-definitivo.csv', 'summary-classic-grasp-RN.csv')
process_csv('Projetos/FJS/classic-SA-definitivo.csv', 'summary-classic-SA.csv')
process_csv('Projetos/FJS/classic-Tayebi.csv', 'summary-classic-Tayebi.csv')
