GlobalOutFile := "../f4(4)_a0a1a3a4_两长两短各一黑一白_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a1,a3,a4]，左侧两长根与右侧两短根各一黑一白";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a1,a3,a4]，左侧两长根与右侧两短根各一黑一白",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a1","a3","a4"],
            target_color_patterns := [
                [1,0,1,0],
                [1,0,0,1],
                [0,1,1,0],
                [0,1,0,1]
            ]
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
