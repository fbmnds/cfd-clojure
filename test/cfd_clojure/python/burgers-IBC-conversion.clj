;; initiallly

;; -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 12.5663706143592)*exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1)))) + 4



;; outer +

(+
;; -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 12.5663706143592)*exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1))))
 4.)



;; '(-8*t + 2*x - 12.5663706143592)' to ' 2.*(-4*t + x - 6.28318530717959)'

(+
;; -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) -  2.*(-4*t + x - 6.28318530717959)*exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1))))
 4.)



;; extract (-4*t + x - 6.28318530717959)

(let [f1 (+ (* -4 t) x -6.28318530717959)]
(+
;; -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) -  2.*f1*exp(-f1**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-f1**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1))))
 4.)
)



;; extract '(-(-4*t + x)**2/f2)'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)]
(+
;; -2*nu*(-(-8*t + 2*x)*exp(f3)/f2 -  2.*f1*exp(-f1**2/f2)/f2)/(exp(-f1**2/f2) + exp(f3))
 4.)
)



;; 'exp(f3)' to '(exp f3)'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)]
(+
;; -2*nu*(-(-8*t + 2*x)*(exp f3)/f2 -  2.*f1*exp(-f1**2/f2)/f2)/(exp(-f1**2/f2) + (exp f3))
 4.)
)



;; extract '-f1**2/f2'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
(+
;; -2*nu*(-(-8*t + 2*x)*(exp f3)/f2 -  2.*f1*exp(f4)/f2)/(exp(f4) + (exp f3))
 4.)
)



;; 'exp(f4)' to '(exp f4)'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
(+
;; -2*nu*(-(-8*t + 2*x)*(exp f3)/f2 -  2.*f1*(exp f4)/f2)/((exp f4) + (exp f3))
 4.)
)


;; reorder outer /, '((exp f4) + (exp f3))'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
(+

 (/ -2*nu*(-(-8*t + 2*x)*(exp f3)/f2 -  2.*f1*(exp f4)/f2)
    (+ (exp f4) (exp f3)))
 4.)
)



;; reorder outer - in /

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
(+

 (/ -2*nu*(-  -(-8*t + 2*x)*(exp f3)/f2
              2.*f1*(exp f4)/f2)
    (+ (exp f4) (exp f3)))
 4.)
)



;; factor out /f2

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
(+

 (/ -2*nu*(-  -(-8*t + 2*x)*(exp f3)
              2.*f1*(exp f4))/f2
    (+ (exp f4) (exp f3)))
 4.)
)



;; reorder /f2

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
(+

 (/ (/ -2*nu*(-  -(-8*t + 2*x)*(exp f3) 2.*f1*(exp f4))
       f2)
    (+ (exp f4) (exp f3)))
 4.)
)




;; reorder -2*nu*(- ...)

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)]
  (+

   (/ (/ (* -2. nu (-  -(-8*t + 2*x)*(exp f3) 2.*f1*(exp f4)))
         f2)
      (+ (exp f4) (exp f3)))
   4.)
  )




;; extract '-(-8*t + 2*x)*(exp f3)'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)
      f5 (* -1. (+ (* -8. t) (* 2. x)) (exp f3))]
  (+
   (/ (/ (* -2. nu (- f5 2.*f1*(exp f4)))
         f2)
      (+ (exp f4) (exp f3)))
   4.)
  )



;; reorder '2.*f1*(exp f4)'

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3 (-(-4*t + x)**2/f2)
      f4 (/ (* -1. f1 f1) f2)
      f5 (* -1. (+ (* -8. t) (* 2. x)) (exp f3))]
  (+
   (/ (/ (* -2. nu (- f5 (* 2. f1 (exp f4))))
         f2)
      (+ (exp f4) (exp f3)))
   4.)
  )


;; reorder f4

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3a (+ (* -4. t) x)
      f3 (-f4**2/f2)
      f4 (/ (* -1. f1 f1) f2)
      f5 (* -1. (+ (* -8. t) (* 2. x)) (exp f3))]
  (+
   (/ (/ (* -2. nu (- f5 (* 2. f1 (exp f4))))
         f2)
      (+ (exp f4) (exp f3)))
   4.)
  )



;; reorder f3

(let [f1 (+ (* -4 t) x -6.28318530717959)
      f2 (* 4. nu (+ t 1.))
      f3a (+ (* -4. t) x)
      f3 (/ (* -1. f4 f4) f2)
      f4 (/ (* -1. f1 f1) f2)
      f5 (* -1. (+ (* -8. t) (* 2. x)) (exp f3))]
  (+
   (/ (/ (* -2. nu (- f5 (* 2. f1 (exp f4))))
         f2)
      (+ (exp f4) (exp f3)))
   4.)
  )
