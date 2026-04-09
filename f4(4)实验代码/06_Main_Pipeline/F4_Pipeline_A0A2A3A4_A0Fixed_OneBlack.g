GlobalOutFile := "../f4(4)_a0a2a3a4_a0固定黑且a2a3a4任一单黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a2,a3,a4]，其中a0固定涂黑，且a2/a3/a4中任一位置涂黑";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a2,a3,a4]，a0固定涂黑，且a2/a3/a4中任意一个涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a2","a3","a4"],
            target_color_patterns := [
                [1,1,0,0],
                [1,0,1,0],
                [1,0,0,1]
            ]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
