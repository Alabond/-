GlobalOutFile := "../f4(4)_a0a1a2a3_任一单黑_实验运行结果.txt";;
GlobalIntroText := "研究对象: F4(4) 实验 -> 选取[a0,a1,a2,a3] 且恰有一个黑点";;

GlobalExperiments := [
    rec(
        desc := "F4(4) 实验: 选取[a0,a1,a2,a3]，四个根中任意一个涂黑",
        roots := rec(
            mode := "SCAN",
            selected_pos_labels := ["a0","a1","a2","a3"],
            required_black_count := 1
        ),
        colors := fail
    )
];;

Read("Main_Pipeline.g");;
