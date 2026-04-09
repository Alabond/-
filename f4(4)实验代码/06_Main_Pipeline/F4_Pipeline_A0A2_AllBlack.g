GlobalOutFile := "../f4(4)_a0a2_全黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a2]，并全部涂黑";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a2]，全部涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a2"],
            target_color_patterns := [
                [1,1]
            ]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
