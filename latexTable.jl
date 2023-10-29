using TexTables


function dm1_latexTable(instances, gl_times, gl_result, ds_times, ds_result)

    t1 = TableCol("time", instances, gl_times)
    t2 = TableCol("Glouton", instances, gl_result)
    t3 = TableCol("time", instances, ds_times)
    t4 = TableCol("Simple Descent", instances, ds_result)

    table = hcat(t1, t2, t3, t4)
    return table
end