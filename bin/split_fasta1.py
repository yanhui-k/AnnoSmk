def split_fasta1(fasta, dir = "R", size = 4000000):
    #file_list = []
    z = 1
    with open(fasta, 'r') as f:
        data = f.read().split('>')
        if os.path.exists(dir) == "True":
            #mv(dir, dir + str(time.time()).replace('.', ''))
            #mkdir(dir)
            pass
        else:
            mkdir(dir)
        z = 1
        file = f"{z}.fa"
        new_file_dir = os.path.join(dir,file)
        open(new_file_dir, "w")
        for i, j in enumerate(data[1:], start=1):
            (header, content) = j.split('\n', 1)
            #new_file_name = f'{header}.fasta'
            #new_file_dir = os.path.join(dir,new_file_name)
            #file_list.append(new_file_name)
            if os.path.getsize(new_file_dir) < size:
                with open(new_file_dir, 'a') as f:
                    f.write('>' + j)
            else:
                z += 1
                file = f"{z}.fa"
                new_file_dir = os.path.join(dir, file)
                with open(new_file_dir,"w") as f:
                    f.write(">" + j)

split_fasta1(snakemake.input[0],snakemake.output[0])