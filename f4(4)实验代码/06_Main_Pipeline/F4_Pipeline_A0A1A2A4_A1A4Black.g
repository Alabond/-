GlobalOutFile := "../f4(4)_a0a1a2a4_a1a4涂黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a1,a2,a4]，固定a1,a4涂黑";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a1,a2,a4] 且 a1,a4 涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a1","a2","a4"],
            target_black_labels := ["a1","a4"]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
