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

(defn burgers-u0 [t x nu]
  (let [f1 (+ (* -4 t) x -6.28318530717959)
        f2 (* 4. nu (+ t 1.))
        f3a (+ (* -4. t) x)
        f3 (/ (* -1. f3a f3a) f2)
        f4 (/ (* -1. f1 f1) f2)
        f5 (* -1. (+ (* -8. t) (* 2. x)) (exp f3))]
    (+ (/ (/ (* -2. nu (- f5 (* 2. f1 (exp f4))))
             f2)
          (+ (exp f4) (exp f3)))
       4.)))
