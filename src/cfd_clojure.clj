(ns cfd-clojure
  (:require [clojure.math.numeric-tower :as math]
            [clojure.core.match :as m]
            [clatrix.core :as clx])
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


(defn discretize-2D [f m u]
  (iterate (partial f m) u))


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


;; Step 5

(defn linear-convection-2D [m un]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        A (sel un :except-rows 0 :except-cols 0)
        B (sel un :except-rows upper_y :except-cols 0)
        C (sel un :except-rows 0 :except-cols upper_x)
        k (* -1. (:c m) (:dt m))
        kx (/ k (:dx m))
        ky (/ k (:dy m))
        u_core (sel (plus (mult (+ 1. kx ky) A)
                          (mult (* -1. kx) B)
                          (mult (* -1. ky) C))
                    :except-rows (dec upper_y)
                    :except-cols (dec upper_x))
        v (matrix 1. (:ny m) (:nx m))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set v y x (sel u_core (dec y) (dec x))))
    v))


;; Step 6:

(defn convection-2D [m [un vn]]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        Au (sel un :except-rows 0 :except-cols 0)
        Bu (sel un :except-rows upper_y :except-cols 0)
        Cu (sel un :except-rows 0 :except-cols upper_x)
        Av (sel vn :except-rows 0 :except-cols 0)
        Bv (sel vn :except-rows upper_y :except-cols 0)
        Cv (sel vn :except-rows 0 :except-cols upper_x)
        k (* -1. (:c m) (:dt m))
        kx (/ k (:dx m))
        ky (/ k (:dy m))
        u_core (sel (plus (mult (+ 1. kx ky) Au)
                          (mult (* -1. kx) Bu)
                          (mult (* -1. ky) Cu))
                    :except-rows (dec upper_y)
                    :except-cols (dec upper_x))
        v_core (sel (plus (mult (+ 1. kx ky) Av)
                          (mult (* -1. kx) Bv)
                          (mult (* -1. ky) Cv))
                    :except-rows (dec upper_y)
                    :except-cols (dec upper_x))
        uu (matrix 1. (:ny m) (:nx m))
        vv (matrix 1. (:ny m) (:nx m))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set uu y x (sel u_core (dec y) (dec x)))
      (clx/set vv y x (sel v_core (dec y) (dec x))))
    [uu vv]))
