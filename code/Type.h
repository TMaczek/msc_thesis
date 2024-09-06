#pragma once

/**
 * @brief 
 * Możliwe typy komórek w programie
 * Niektóre typy opisane w pracy Harlova-Welsha (corner, obstacle, URON) zostały zrealizowane przez funkcje w klasie System.
 * Dzięki temu każda komórka ma jedynie jeden typ (a nie więcej).
 * Zgrupowania typów (np. brzegowe, z cieczą) zrealizowane są w Cell
 */
enum Type{ 
	border = 1, 		// Brzeg lub przeszkoda, domyślnie z poślizgiem
	full = 2,			// Komórka pełna cieczy
	surface = 3,  		// Komórka powierzchniowa cieczy
	emptyc = 4, 		// Komórka pusta (próżnia/powietrze) (emptyc(ell))
	border_noslip = 5, 	// Brzeg lub przeszkoda bez poślizgu
	out = 6				// Komórka wyjściowa
};						// Można dorobić komórkę wejściową (należy dodać jej implementację w System)