
def imprimir_gantt(arquivo_entrada, arquivosaida):
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.font_manager import FontProperties

    df = pd.read_csv(arquivo_entrada, header=None, names=["Task", "Start", "Finish", "Resource", "Operation"] )
    df["Diff"] = df.Finish - df.Start
    
    color = {"Type 1":"0.25", "Type 2":"0.5", "Type 3":"0.75", "Type 0":"1", "Type 4":"blue", "Type 5":"red", "Type 6":"green"}
    color2 = {"Type 1":"0", "Type 2":"0", "Type 3":"0", "Type 0":['none'], "Type 4":"0", "Type 5":"0", "Type 6":"0"}
    
    
    blackpatch	= mpatches.Patch(facecolor='0.25', label='Job 1')
    graypatch	= mpatches.Patch(color='0.5', label='Job 2')
    whitepatch	= mpatches.Patch(color='0.75', label='Job 3')
    bluepatch	= mpatches.Patch(color='blue', label='Job 4')
    #redpatch	= mpatches.Patch(color='red', label='Job 5')
    #greenpatch   = mpatches.Patch(color='green', label='Job 6')
    
    fig,ax=plt.subplots(figsize=(10,2))
    
    
    labels= ['M1', 'M2', 'M3', 'M4', 'M5'] #,'M6','M7']
    
    for i, task in enumerate(df.groupby("Task")):
        labels.append(task[0])
        for r in task[1].groupby("Resource"):
            data = r[1][["Start", "Diff"]]
            ax.broken_barh(data.values, (i-0.4,0.8), facecolor=color[r[0]], edgecolor=color2[r[0]] )

            data2 = r[1][["Start", "Diff","Operation"]]
            #print(data2)
            
            for x1, x2, x3 in data2.values:
                if x3 == -1 or x2 < 5 :
                    continue
                ax.text(x=x1 + x2/2, 
                        y=i,
                        s= int(x3), 
                        ha='center', 
                        va='center',
                        color='white',
                    )

    
    
    
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels) 
    #ax.set_ylim([-0.5, 4.5])
    ax.set_ylim([-0.5, 4.5])
    ax.set_xlabel("Time")
    
    fontP = FontProperties()
    fontP.set_size('small')
    #plt.legend(handles=[blackpatch,graypatch,whitepatch,bluepatch],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.legend(handles=[blackpatch,graypatch,whitepatch,bluepatch],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.1)


    plt.tight_layout()       
    #plt.show()
    plt.savefig(arquivosaida, bbox_inches='tight', dpi = 250)