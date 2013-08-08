(ns cfd-clojure
  (:require [clojure.math.numeric-tower :as math]
            [clojure.core.match :as m])
  (:use (incanter core charts)))

;; provide discretisation linear in t and up to second order in x


(declare burgers-bc-err)
(defn discretize [f m u]
  (let [g0 (fn [u] (partition (:x-steps m) 1 u))
        h0 (fn [v] (map #(f m %) v))
        l (last u)
        g (fn [u] (conj (map #(f m %) (partition (:x-steps m) 1 u)) (first u)))
        h (fn [u] ((comp g #(concat % [l])) u))]
    (m/match [(:x-steps m) (:bc m)]
             [2 nil] (iterate g u)
             [3 nil] (iterate h u)
             [3 burgers-bc-err] (iterate (comp h0 (:bc m) g0) u)
             [3 _] (iterate (comp h0 g0 (:bc m)) u)
             ;; :else throw IllegalArgumentException
             )))



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


;; In the given model, u0 and u101 refer to the same location modulo 2*pi.
;; Hence, the same location is used twice in updating the BC.
;; This is most likely a bug.
;;
(defn burgers-bc-err [g0u]
  (let [[u0 u1 u2] (vec (first g0u))
        [u99 u100 u101] (vec (last g0u))]
    (concat (conj g0u [u101 u0 u1]) [[u100 u101 u0]])))


;; (~ denotes 2*pi periodicy)
;; if un[0] ~ un[nx]
;; then un[1] ~ un[nx + 1]
;; and un[-1] ~ un[nx - 1]
;;
(defn burgers-bc [u]
  (concat (conj u (last (butlast u))) [(second u)]))


(defn burgers-eqn [m [u_ni-1 u_ni u_ni+1]]
  (let [dx (:dx m)
        nu (:nu m)
        dt (:dt m)]
    (+ u_ni
       (* u_ni (/ dt dx) (- u_ni-1 u_ni))
       (* nu (/ dt (* dx dx)) (+ u_ni+1 (* -2. u_ni) u_ni-1)))))










(comment
(defn burgers [f m u]
  (let [un+1 (second u)
        un-0 (first u)
        un-1 (last u)
        un-2 (last (butlast u))
        u_bc (fn [un-0 un-1 un-2] (+ un-1
                                    (* un-1 (/ dt dx) (- un-2 un-1))
                                    (* nu (/ dt (* dx dx)) (+ un-0 (* -2. un-1) un-2))))
        u-1 (u_bc un-0 un-1 un-2)
        u-0 (u_bc un+1 un-0 un-1)
        g (fn [u] (conj (map #(f m %) (partition (:x-steps m) 1 u)) u-0))]
    (g (concat u [u-1]))))
)


;; un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*(un[0]-2*un[-1]+un[-2])
