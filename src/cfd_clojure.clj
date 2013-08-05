(ns cfd-clojure
  (:require [clojure.math.numeric-tower :as math])
  (:use (incanter core charts)))

;; provide discretisation linear in t and up to second order in x

(defn discretize [f m u]
  (let [l (last u)
        g (fn [u] (conj (map #(f m %) (partition (:x-steps m) 1 u)) (first u)))
        h (fn [u] (g (concat u [l])))]
    (if (= 2 (:x-steps m))
      (iterate g u)
      (iterate h u))))

;; Step 1

(defn linear-convection [m [u_ni-1 u_ni]]
  (- u_ni (/ (* (:c m) (:dt m) (- u_ni u_ni-1)) (:dx m))))

;; Step 2

(defn non-linear-convection [m [u_ni-1 u_ni]]
  (- u_ni (/ (* u_ni (:dt m) (- u_ni u_ni-1)) (:dx m))))

;; Step 3

(defn diffusion-1D [m [u_ni-1 u_ni u_ni+1]]
  (+ u_ni (/ (* (:nu m) (:dt m) (+ u_ni+1 u_ni-1 (* -2. u_ni))) (Math/pow (:dx m) 2.))))


;; Step 4

;; -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) - (-8*t + 2*x - 12.5663706143592)*exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 6.28318530717959)**2/(4*nu*(t + 1))) + exp(-(-4*t + x)**2/(4*nu*(t + 1)))) + 4

;;(defn- exp [x]
;;  (math/expr Math/E x))

;;(defn burgers-u0 [t x nu]
;;  (let [f1 (exp -1.*(math/expr (+ (* -4. t) x) 2.)/(* 4. nu (inc t)))
;;        f2 (* 4. nu (inc t))
;;        f3 (math/expr (+ (* -4. t) x -6.28318530717959) 2.)
;;        f4 (* 2. (+ (* -4. t) x -6.28318530717959))
;;        f5 (exp -1.*f3/f2)
;;        f6 (* -1. (+ (* -8. t) (* 2. x)))]
;;    (+ (/ (* -2. nu (- f6*f1/f2 f4*f5/f2) (+ f5 f1)) 4.))))
