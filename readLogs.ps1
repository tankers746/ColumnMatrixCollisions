write-host "nodes, ppn, walltime"
for($i = 3625; $i -le 3774; $i++) {
    $filename = ".\logs\auto.sh.o" + $i
    $logfile = Get-Content $filename

    $wall = $logfile | select-string "Wall" | %{$_.line}
    write-host -NoNewline ($logfile | select-string "Nodes" | %{$_.line}).split("Nodes :").split("PPN")
    if($wall) {
        write-host "," $wall.split("Wall time      :").replace("                  ", "")
    } else {
        write-host ",0"
    }

    
}