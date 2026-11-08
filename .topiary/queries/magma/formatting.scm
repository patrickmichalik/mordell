; Formatting for the expression-based Magma check.
[(string) (comment)] @leaf
[(assignment) (print_statement) (comment)] @allow_blank_line_before

; Separate statements and preserve paragraph breaks.
(program ";" @append_hardline)
";" @prepend_antispace
(comment) @prepend_input_softline @append_hardline
(print_statement "print" @append_space)

; Keep function calls and generator lists compact.
(argument_list "(" @prepend_antispace @append_antispace)
")" @prepend_antispace
(constructor "<" @prepend_antispace @append_antispace)
(constructor ">" @prepend_antispace)
"," @prepend_antispace @append_space
":" @prepend_space @append_space
":=" @prepend_space @append_space

; Space arithmetic operators except exponentiation and unary signs.
(binary_operator
  operator: ["+" "-" "*" "/"] @prepend_space @append_space)
(binary_operator
  operator: "^" @prepend_antispace @append_antispace)
(unary_operator
  operator: ["+" "-"] @append_antispace)
